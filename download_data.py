import os
import ncbi.datasets
import pandas as pd 
from datetime import datetime
from zipfile import ZipFile
from Bio import SeqIO
import gffutils
import shutil
import subprocess
import gzip
import urllib3
import argparse

def build_arg_parser():
    """Handle the parameter sent when executing the script from
 the terminal

    Returns
    -----------
    A parser object with the chosen option and parameters"""

    parser = argparse.ArgumentParser(description="This script download data from UniProt Reference Proteomes, NCBI or Ensembl databases and organizes it in a folder structure compatible with OMArk's Snakemake implementation.")
    parser.add_argument('-d', '--directory', type=str, help="Path to the output directory.", required=True)
    parser.add_argument('-db', '--database', type=str, help='Database to download from', choices=['UniProt', 'NCBI', 'Ensembl', 'All'], required=True)
    parser.add_argument('-t', '--datetime',  type=datetime.fromisoformat, help='Date from last update in ISOformat - YYYY-MM-DD', default=None)
    return parser


def get_all_assemblies_taxon_NCBI(tax_name,api_instance):
    assemblies = []
    token = ''
    while token!=None:
        genome_summary = api_instance.assembly_descriptors_by_taxon(
            taxon=tax_name,
            filters_reference_only = True,
            filters_has_annotation = True,
            page_token=token,
            page_size=1000)
        token = genome_summary.next_page_token
        assemblies+= map(lambda d: d.assembly, genome_summary.assemblies)
    return assemblies
    

def download_all_assemblies_NCBI(assemblies, work_path,api_instance, date=None, step=50, restart=None):
    for i in range(0, len(assemblies), step):    
        subset_assemblies = assemblies[i:i+step]
        all_accessions = []
        if restart and i < restart:
            continue
        for assembly in subset_assemblies:

            date_str = assembly.annotation_metadata.release_date
            date_annotation = pd.to_datetime(date_str)

            if len(date_str.split('-'))==3:
                date_annotation = datetime.strptime(date_str,'%Y-%m-%d')
            elif len(date_str.split('-'))==2:
                date_annotation = datetime.strptime(date_str,'%Y-%m')
            else:
                date_annotation = datetime.strptime(date_str,'%Y')
            if not date or date_annotation>date:
                all_accessions.append(assembly.assembly_accession)
        exclude_sequence = True
        include_annotation_type = ['PROT_FASTA', 'GENOME_GFF']
        
        api_response = api_instance.download_assembly_package(
            all_accessions,
            exclude_sequence=exclude_sequence,
            include_annotation_type=include_annotation_type,
            # Because we are streaming back the results to disk, 
            # we should defer reading/decoding the response
            _preload_content=False
        )

        with open(work_path+f'/NCBI_ref_proteomes_{i}.zip', 'wb') as f:
            f.write(api_response.data)
        with ZipFile(work_path+f'/NCBI_ref_proteomes_{i}.zip', 'r') as zObject:
            zObject.extractall(path=work_path)


def extract_data_NCBI( data_dir, fasta_dir, splice_dir, resume=True):
	for folder in os.listdir(data_dir):

	    gene_to_prot= {}
	    prot_to_gene = {}
	    folder_path = data_dir+'/'+folder
	    if os.path.isdir(folder_path):
	        gff = folder_path+'/genomic.gff'
	        if not os.path.isfile(gff):
	            print('Missing file - '+gff)
	        if resume and os.path.isfile(splice_dir+"/"+folder+'.splice') and os.path.isfile(fasta_dir+'/'+folder+'.fa'):
	        	continue
	        db = gffutils.create_db(gff, ':memory:', merge_strategy="create_unique", keep_order=True)
	        for t in db.features_of_type('gene', order_by='start'):
	            gene = t.id
	            gene_list = []
	            ordered_child = list(db.children(t, featuretype='CDS', order_by='start'))
	            for child in ordered_child:
	                type_attribute = ['protein_id', 'Name']
	                for att_type in type_attribute:
	                    protein = child.attributes.get(att_type, [None])[0]
	                    if protein:

	                        break
	                if not protein:
	                    print('warning')
	                    print(child)
	                    continue
	                corr_gene = prot_to_gene.get(protein, None)
	                if corr_gene and corr_gene!=gene:
	                    gene = corr_gene
	                    for other_prot in gene_list:
	                        prot_to_gene[other_prot] = gene
	                    gene_list = gene_to_prot[gene]+gene_list
	                else:
	                    prot_to_gene[protein] = gene
	                if protein not in gene_list:
	                    gene_list.append(protein)
	            if len(gene_list)!=0:
	                gene_to_prot[gene] = gene_list
	        shutil.copyfile(folder_path+'/protein.faa', fasta_dir+'/'+folder+'.fa')
	    write_splice_file(splice_dir+"/"+folder+'.splice','w') 

def write_splice_file(splice_data, splice_file):
    with open(splice_file,'w') as handle_output:
        for val in splice_data.values():
            handle_output.write(";".join(val)+'\n')

def prepare_data_ensembl(fasta_file, splice_dir):
		print(fasta_file)
		splice_value = extract_splice_data_ensembl(fasta_file)
		write_splice_file(splice_value, os.path.join(splice_dir, ".".join(os.path.basename(fasta_file).split('.')[:-1]+['splice'])))


def prepare_all_files_ensembl(fasta_dir, splice_dir):
        for fasta_file in os.listdir(fasta_dir):
            #if not os.path.isfile(os.path.join(fasta_dir, folder+'.fa')):
            prepare_data_ensembl(os.path.join(fasta_dir, fasta_file), splice_dir)


def store_metadata_NCBI(assemblies, work_path):
	data = []
	for assembly in assemblies:

		accession = assembly.assembly_accession
		taxid = assembly.org.tax_id
		date_annotation = assembly.annotation_metadata.release_date
		data.append([accession, taxid, date_annotation])
	with open(work_path+"/GC_to_taxid.csv",'w') as of:
		for assembly in data :
			of.write('\t'.join(assembly)+'\n')

def make_folder_struct(folder):
	if not os.path.isdir(folder):
		os.mkdir(folder)
	for extension in ['Source', 'Splice', 'tmp']:
		if not os.path.isdir(os.path.join(folder, extension)):
			os.mkdir(os.path.join(folder, extension))

def extract_splice_data_ensembl(fasta_file):
    all_splice = {}

    with open(fasta_file, "r") as handle:
        
        fasta_read = SeqIO.parse(handle, "fasta")
        for seq_record in SeqIO.parse(handle, "fasta"):
            string_to_split = seq_record.description
            split = string_to_split.split('gene:') [1]
            gene_id = split.split(' ')[0]
            all_splice[gene_id] = all_splice.get(gene_id, [])+[seq_record.id]
    return all_splice    

def download_all_proteomes_Ensembl(fasta_folder):

	all_ens_ftp = [f"wget -r -nd  -P {fasta_folder} -A '*.pep.all.fa.gz' ftp://ftp.ensembl.org/pub/current_fasta/"	,
	f"wget -r -nd -P {fasta_folder} -A '*.pep.all.fa.gz' ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/current/plants/fasta/",
	f"wget -r -nd -P {fasta_folder} -A '*.pep.all.fa.gz' ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/",
	f"wget -r -nd  -P {fasta_folder} -A '*.pep.all.fa.gz' ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/current/protists/fasta/",
	f"wget -r -nd -P {fasta_folder} -A '*.pep.all.fa.gz' ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/current/fungi/fasta/"]
	for ftp_repo in all_ens_ftp:
		process = subprocess.Popen(
	        ftp_repo,
	        stdout = subprocess.PIPE,
	        stderr = subprocess.PIPE,
	        text= True,
	        shell = True
	    )
		std_out, std_err = process.communicate()
	for file in os.listdir(fasta_folder):
		full_gz_file = os.path.join(fasta_folder, file)
		if file[-2:]=='gz':
			with gzip.open(full_gz_file, 'rb') as f_in:
			    with open(full_gz_file[:-3], 'wb') as f_out:
			        shutil.copyfileobj(f_in, f_out)
			os.remove(full_gz_file)
		elif file[-4:-1]=='gz.':
			os.remove(full_gz_file)
def download_metadata_ensembl(folder):
	files_url = ['ftp.ensembl.org/pub/current/species_EnsemblVertebrates.txt', 'ftp.ensemblgenomes.org/pub/metazoa/current/species_EnsemblMetazoa.txt' , 'ftp.ensemblgenomes.org/pub/plants/current/species_EnsemblPlants.txt' , 
	'ftp.ensemblgenomes.org/pub/protists/current/species_EnsemblProtists.txt', 'ftp.ensemblgenomes.org/pub/fungi/current/species_EnsemblFungi.txt']
	full_info = []
	http = urllib3.PoolManager()

	for url in files_url:
		resp = http.request("GET", url)
		text_answer = resp.data.decode()
		for line in text_answer.split('\n'):
			if line!='' and line[0]!='#':
				full_info.append(line)
	outfile = os.path.join(folder, 'Info.txt')
	with open(outfile, 'w') as of:
		of.write('\n'.join(full_info))


def download_all_proteomes_UP(fasta_folder):
	command = f"wget -r -nd  -P {fasta_folder} -A '*.fasta.gz' ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/"
	process = subprocess.Popen(
		       command,
		       stdout = subprocess.PIPE,
		       stderr = subprocess.PIPE,
		       text= True,
		       shell = True
		    )
	std_out, std_err = process.communicate()
	for file in os.listdir(fasta_folder):
		if file[-2:]=='gz':
			full_gz_file = os.path.join(fasta_folder, file)

			if 'DNA' not in file and 'additional' not in file:
				with gzip.open(full_gz_file, 'rb') as f_in:
					f_out_name = full_gz_file[:-3]
					f_out_name = f_out_name.replace('fasta', 'fa')
					with open(f_out_name, 'wb') as f_out:
						shutil.copyfileobj(f_in, f_out)
			os.remove(full_gz_file)



if __name__ =='__main__':

	parser = build_arg_parser()
	arg = parser.parse_args()
	directory = arg.directory
	source = arg.database
	date = arg.datetime
	make_folder_struct(directory)
	tmp_folder = os.path.join(directory,'tmp')
	fasta_folder = os.path.join(directory,'Source')
	splice_folder= os.path.join(directory,'Splice')

	if source =='NCBI':
		api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())
		assemblies = get_all_assemblies_taxon_NCBI('Eukaryota',api_instance)
		download_all_assemblies_NCBI(assemblies, tmp_folder,api_instance)
		store_metadata_NCBI(assemblies, tmp_folder)
		extract_data_NCBI(os.path.join(tmp_folder,'ncbi_dataset', 'data'), fasta_folder, splice_folder, date=date)
	elif source == 'Ensembl':
		download_all_proteomes_Ensembl(fasta_folder)
		prepare_all_files_ensembl(fasta_folder, splice_folder)
		download_metadata_ensembl(directory)
	elif source == 'UniProt':
		download_all_proteomes_UP(fasta_folder)
