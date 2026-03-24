#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 19:28:33 2022

@author: omidard
@updater: wwzyeo

Generate draft genome-scale metabolic models (GEMs) from a reference model and
target genomes. The pipeline builds orthology matrices using BLAST, prunes a
reference model per strain, and writes per-strain draft GEMs as JSON.

Sample usage:
    python GEMgenerator.py \
      --reactome LBReactome.json \
      --genes-fna data/sequences/genes.fna \
      --genes-faa data/sequences/genes.faa \
      --genbank  data/genomes/ \
      --out      results/gems_draft \
      --threads  96

Key inputs (defaults shown):
    --reactome   LBReactome.json            (SBML/XML or JSON reference model)
    --genes-fna  reactome_genes/reactome_nucleotide.fa  (reference nucleotide FASTA)
    --genes-faa  reactome_genes/reactome_prot.fa        (reference protein FASTA)
    --genbank    target_genome_dir         (directory of target .gbk/.gbff files)
    --out        .                         (base output directory)
    --threads    16                        (BLAST threads)
"""
import argparse
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.delete import delete_model_genes, remove_genes
import os
from os.path import join
import pandas as pd
from glob import glob
from Bio import Entrez, SeqIO
from os.path import join
from os.path import isfile, join
from os import listdir
import shutil
import numpy as np



def get_file_names (directory1,directory2):
    onlyfiles = [f for f in listdir(directory1) if isfile(join(directory1, f)) and (f.endswith(".gbk") or f.endswith(".gbff"))]
    gl8 = []
    for i in onlyfiles:
        if i.endswith(".gbff"):
            gl1 = i.replace(".gbff", "")
        else:
            gl1 = i.replace(".gbk", "")
        gl8.append(gl1)
    gl9=pd.DataFrame({'strain':gl8,'NCBI ID':gl8, 'Pathotype':range(len(onlyfiles))})
    path = directory2+'/StrainInformation.xlsx'
    StrainInformation = gl9.to_excel(path)


def dl_genome(id, folder='genomes'): # be sure get CORRECT ID
    files=glob('%s/*.gbk'%folder)
    out_file = '%s/%s.gbk'%(folder, id)

    if out_file in files:
        print (out_file, 'already downloaded')
        return
    else:
        print ('downloading %s from NCBI'%id)
        
    from Bio import Entrez
    Entrez.email = "omidard@biosustain.dtu.dk"     #Insert email here for NCBI
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
    fout = open(out_file,'w')
    fout.write(handle.read())
    fout.close()


def get_strain_info(directory1):
    files = glob('%s/*.gbk'%directory1) + glob('%s/*.gbff'%directory1)
    strain_info = []
    
    for file in files:
        handle = open(file)
        record = SeqIO.read(handle, "genbank")
        for f in record.features:
            if f.type=='source':
                info = {}
                info['file'] = file
                info['id'] = file.split('\\')[-1].split('.')[0]
                for q in f.qualifiers.keys():
                    info[q] = '|'.join(f.qualifiers[q])
                strain_info.append(info)
    return pd.DataFrame(strain_info)


def resolve_genbank_path(id, folder):
    gbk = '%s/%s.gbk' % (folder, id)
    gbff = '%s/%s.gbff' % (folder, id)
    if os.path.isfile(gbk):
        return gbk
    if os.path.isfile(gbff):
        return gbff
    return gbk


def parse_genome(id, type='prot', in_folder='genomes', out_folder='prots', overwrite=1):

    in_file = resolve_genbank_path(id, in_folder)
    out_file='%s/%s.fa'%(out_folder, id)
    files =glob('%s/*.fa'%out_folder)
    
    if out_file in files and overwrite==0:
        print (out_file, 'already parsed')
        return
    else:
        print ('parsing %s'%id)
    
    handle = open(in_file)
    
    fout = open(out_file,'w')
    x = 0
    
    records = SeqIO.parse(handle, "genbank")
    for record in records:
        for f in record.features:
            if f.type=='CDS':
                seq=f.extract(record.seq)
                
                if type=='nucl':
                    seq=str(seq)
                else:
                    seq=str(seq.translate())
                    
                if 'locus_tag' in f.qualifiers.keys():
                    locus = f.qualifiers['locus_tag'][0]
                elif 'gene' in f.qualifiers.keys():
                    locus = f.qualifiers['gene'][0]
                else:
                    locus = 'gene_%i'%x
                    x+=1
                fout.write('>%s\n%s\n'%(locus, seq))
    fout.close()


def make_blast_db(id,folder='prots',db_type='prot'):
    import os
    
    out_file ='%s/%s.fa.pin'%(folder, id)
    files =glob('%s/*.fa.pin'%folder)
    
    if out_file in files:
        print (id, 'already has a blast db')
        return
    if db_type=='nucl':
        ext='fna'
    else:
        ext='fa'

    cmd_line="makeblastdb -in %s/%s.%s -dbtype %s" %(folder, id, ext, db_type)
    
    print (' making blast db with following command line...')
    print (cmd_line)
    os.system(cmd_line)

def run_blastp(seq,db,in_folder='prots_dir', out_folder='bbh_dir', out=None,outfmt=6,evalue=0.001,threads=64):
    import os
    if out==None:
        out='%s/%s_vs_%s.txt'%(out_folder, seq, db)
        print(out)
    
    files =glob('%s/*.txt'%out_folder)
    if out in files:
        print (seq, 'already blasted')
        return
    
    print ('blasting %s vs %s'%(seq, db))
    
    db = '%s/%s.fa'%(in_folder, db)
    seq = '%s/%s.fa'%(in_folder, seq)
    cmd_line='blastp -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' %(db, seq, out, evalue, outfmt, threads)
    
    print ('running blastp with following command line...')
    print (cmd_line)
    os.system(cmd_line)
    return out


def get_gene_lens(query, in_folder='prots'):

    file = '%s/%s.fa'%(in_folder, query)
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = []
    
    for record in records:
        out.append({'gene':record.name, 'gene_length':len(record.seq)})
    
    out = pd.DataFrame(out)
    return out


def get_bbh(query, subject, in_folder='bbh', prots_dir='prots_dir', threads=64):    
    
    #Utilize the defined protein BLAST function
    run_blastp(query, subject, in_folder=prots_dir, out_folder=in_folder, threads=threads)
    run_blastp(subject, query, in_folder=prots_dir, out_folder=in_folder, threads=threads)
    
    query_lengths = get_gene_lens(query, in_folder=prots_dir)
    subject_lengths = get_gene_lens(subject, in_folder=prots_dir)
    
    #Define the output file of this BLAST
    out_file = '%s/%s_vs_%s_parsed.csv'%(in_folder,query, subject)
    files=glob('%s/*_parsed.csv'%in_folder)
    
    #Combine the results of the protein BLAST into a dataframe
    print ('parsing BBHs for', query, subject)
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    bbh=pd.read_csv('%s/%s_vs_%s.txt'%(in_folder,query, subject), sep='\t', names=cols)
    bbh = pd.merge(bbh, query_lengths) 
    bbh['COV'] = bbh['alnLength']/bbh['gene_length']
    
    bbh2=pd.read_csv('%s/%s_vs_%s.txt'%(in_folder,subject, query), sep='\t', names=cols)
    bbh2 = pd.merge(bbh2, subject_lengths) 
    bbh2['COV'] = bbh2['alnLength']/bbh2['gene_length']
    out = pd.DataFrame()
    
    # Filter the genes based on coverage
    bbh = bbh[bbh.COV>=0.30]
    bbh2 = bbh2[bbh2.COV>=0.30]
    
    # Vectorized BBH selection: replace per-gene loop with groupby+idxmax for speed.
    if not bbh.empty and not bbh2.empty:
        bbh_best = bbh.loc[bbh.groupby("gene")["PID"].idxmax()].copy()
        bbh2_best = bbh2.loc[bbh2.groupby("gene")["PID"].idxmax()][["gene", "subject"]].copy()
        bbh2_map = bbh2_best.set_index("gene")["subject"]
        bbh_best["BBH"] = np.where(bbh_best["subject"].map(bbh2_map) == bbh_best["gene"], "<=>", "->")
        out = bbh_best
    
    #Save the final file to a designated CSV file
    out.to_csv(out_file)


def gbk2fasta(gbk_filename):
    faa_filename = '.'.join(gbk_filename.split('.')[:-1])+'.fna'
    input_handle  = open(gbk_filename, "r")
    output_handle = open(faa_filename, "w")

    for seq_record in SeqIO.parse(input_handle, "genbank") :
        print ("Converting GenBank record %s" % seq_record.id)
        output_handle.write(">%s %s\n%s\n" % (
               seq_record.id,
               seq_record.description,
               seq_record.seq))

    output_handle.close()
    input_handle.close()


def run_blastn(seq, db, nucls_dir='nucls_dir', target_genome_dir='target_genome_dir', outfmt=6, evalue=0.001, threads=64):
    import os
    out = nucls_dir+'/'+seq+'_vs_'+db+'.txt'
    seq = nucls_dir+'/'+seq+'.fa'
    db = target_genome_dir+'/'+db+'.fna'
    
    cmd_line='blastn -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' \
    %(db, seq, out, evalue, outfmt, threads)
    
    print ('running blastn with following command line...')
    print (cmd_line)
    os.system(cmd_line)
    return out



def parse_nucl_blast(infile):
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    data = pd.read_csv(infile, sep='\t', names=cols)
    data = data[(data['PID']>50) & (data['alnLength']>0.5*data['queryEnd'])]
    data2=data.groupby('gene').first()
    return data2.reset_index()


def extract_seq(g, contig, start, end):
    from Bio import SeqIO
    handle = open(g)
    records = SeqIO.parse(handle, "fasta")
    
    seq = None
    for record in records:
        if record.name==contig:
            if end>start:
                section = record[start:end]
            else:
                section = record[end-1:start+1].reverse_complement()
                
            seq = str(section.seq)
    if seq is None:
        raise ValueError(f"Contig '{contig}' not found in {g}")
    return seq
def load_reference_model(path):
    if path.endswith(".json"):
        return cobra.io.load_json_model(path)
    return cobra.io.read_sbml_model(path)


def main():
    parser = argparse.ArgumentParser(description="Generate draft GEMs from genomes and a reference model.")
    parser.add_argument("--reactome", default="LBReactome.json", help="Reference model (SBML/XML or JSON).")
    parser.add_argument("--genes-fna", default="reactome_genes/reactome_nucleotide.fa", help="Reference nucleotide sequences (FASTA).")
    parser.add_argument("--genes-faa", default="reactome_genes/reactome_prot.fa", help="Reference protein sequences (FASTA).")
    parser.add_argument("--genbank", default="target_genome_dir", help="Directory of target GenBank files (.gbk/.gbff).")
    parser.add_argument("--out", default=".", help="Base output directory for generated folders.")
    parser.add_argument("--threads", type=int, default=16, help="Threads for BLAST.")
    args = parser.parse_args()

    base_out = args.out
    os.makedirs(base_out, exist_ok=True)

    gap_inf_dir = join(base_out, "gap_inf_dir")
    reference_genome_dir = join(base_out, "reference_genome_dir")
    target_genome_dir = args.genbank
    prots_dir = join(base_out, "prots_dir")
    nucls_dir = join(base_out, "nucls_dir")
    bbh_dir = join(base_out, "bbh_dir")
    present_absence_dir = join(base_out, "present_absence_dir")
    initial_models_dir = join(base_out, "initial_models_dir")
    output_models_dir = join(base_out, "output_models_dir")
    gapfilled_models_dir = join(base_out, "gapfilled_models_dir")
    blast_exe_dir = join(base_out, "blast_exe_dir")
    temp_files = join(base_out, "temp_files")
    ref_model_dir = join(base_out, "ref_model_dir")

    directories = [
        gap_inf_dir,
        reference_genome_dir,
        prots_dir,
        nucls_dir,
        bbh_dir,
        present_absence_dir,
        initial_models_dir,
        output_models_dir,
        gapfilled_models_dir,
        blast_exe_dir,
        temp_files,
        ref_model_dir,
    ]
    for d in directories:
        os.makedirs(d, exist_ok=True)

    if args.genes_faa:
        shutil.copyfile(args.genes_faa, join(prots_dir, "reactome.fa"))
    if args.genes_fna:
        shutil.copyfile(args.genes_fna, join(nucls_dir, "reactome.fa"))
    if not args.genes_faa and not os.path.isfile(join(prots_dir, "reactome.fa")):
        raise FileNotFoundError(
            "Missing reference protein FASTA. Provide --genes-faa or place reactome.fa in prots_dir."
        )
    if not args.genes_fna and not os.path.isfile(join(nucls_dir, "reactome.fa")):
        raise FileNotFoundError(
            "Missing reference nucleotide FASTA. Provide --genes-fna or place reactome.fa in nucls_dir."
        )

    get_file_names(target_genome_dir, temp_files)
    StrainsOfInterest = pd.read_excel(join(temp_files, "StrainInformation.xlsx"))
    print(StrainsOfInterest)
    referenceStrainID = "reactome"
    targetStrainIDs = list(StrainsOfInterest["NCBI ID"])

    files = glob("%s/*.gbk" % target_genome_dir) + glob("%s/*.gbff" % target_genome_dir)
    strain_info = []
    for file in files:
        handle = open(file)
        record = list(SeqIO.parse(handle, "genbank"))
        for i in record:
            for f in i.features:
                if f.type == "source":
                    info = {}
                    info["file"] = file.replace(target_genome_dir, "")
                    info["id"] = file.replace(target_genome_dir, "")
                    info["genome_size"] = len(i.seq) / 1000000
                    for q in f.qualifiers.keys():
                        info[q] = "|".join(f.qualifiers[q])
                        strain_info.append(info)
    sinf = pd.DataFrame(strain_info)
    sinf2 = sinf.drop_duplicates(subset="id", keep="first", inplace=False, ignore_index=False)
    print(sinf2)

    for strain in targetStrainIDs:
        parse_genome(strain, type="prot", in_folder=target_genome_dir, out_folder=prots_dir)
        parse_genome(strain, type="nucl", in_folder=target_genome_dir, out_folder=nucls_dir)

    for strain in targetStrainIDs:
        make_blast_db(strain, folder=prots_dir, db_type="prot")
    make_blast_db(referenceStrainID, folder=prots_dir, db_type="prot")

    for strain in targetStrainIDs:
        get_bbh(referenceStrainID, strain, in_folder=bbh_dir, prots_dir=prots_dir, threads=args.threads)

    blast_files = glob("%s/*_parsed.csv" % bbh_dir)
    for blast in blast_files:
        bbh = pd.read_csv(blast)
        print(blast, bbh.shape)

    if not os.path.isfile(args.reactome):
        reactome_alt = "LBReactome.xml"
        if os.path.isfile(reactome_alt):
            print(f"Reference model not found at {args.reactome}; using {reactome_alt} instead.")
            args.reactome = reactome_alt
        else:
            raise FileNotFoundError(f"Reference model not found: {args.reactome}")
    model = load_reference_model(args.reactome)
    listGeneIDs = []
    for gene in model.genes:
        listGeneIDs.append(gene.id)
    ref_gene_count = len(model.genes)
    ref_rxn_count = len(model.reactions)

    ortho_matrix = pd.DataFrame(index=listGeneIDs, columns=targetStrainIDs)
    geneIDs_matrix = pd.DataFrame(index=listGeneIDs, columns=targetStrainIDs)
    print(len(listGeneIDs))

    for blast in blast_files:
        bbh = pd.read_csv(blast)
        listIDs = []
        listPID = []
        for r, row in ortho_matrix.iterrows():
            try:
                currentOrtholog = bbh[bbh["gene"] == r].reset_index()
                listIDs.append(currentOrtholog.iloc[0]["subject"])
                listPID.append(currentOrtholog.iloc[0]["PID"])
            except:
                listIDs.append("None")
                listPID.append(0)
        for col in ortho_matrix.columns:
            if col in blast:
                ortho_matrix[col] = listPID
                geneIDs_matrix[col] = listIDs
    print(sorted(ortho_matrix))

    for column in ortho_matrix:
        ortho_matrix.loc[ortho_matrix[column] <= 50.0, column] = 0
        ortho_matrix.loc[ortho_matrix[column] > 50.0, column] = 1

    for strain in targetStrainIDs:
        gbk2fasta(resolve_genbank_path(strain, target_genome_dir))

    for strain in targetStrainIDs:
        make_blast_db(strain, folder=target_genome_dir, db_type="nucl")

    genome_blast_res = []
    for strain in targetStrainIDs:
        res = run_blastn(referenceStrainID, strain, nucls_dir=nucls_dir, target_genome_dir=target_genome_dir, threads=args.threads)
        genome_blast_res.append(res)

    na_rows = []
    for file in genome_blast_res:
        genes = parse_nucl_blast(file)
        name = ".".join(file.split("_")[-1].split(".")[:-1])
        na_rows.append(genes[["gene", "subject", "PID"]])
    na_matrix = pd.concat(na_rows, ignore_index=True)
    na_matrix = pd.pivot_table(na_matrix, index="gene", columns="subject", values="PID")

    ortho_matrix_w_unannotated = ortho_matrix.copy()
    geneIDs_matrix_w_unannotated = geneIDs_matrix.copy()

    nonModelGenes = []
    for g in na_matrix.index:
        if g not in listGeneIDs:
            nonModelGenes.append(g)

    na_model_genes = na_matrix.drop(nonModelGenes)

    pseudogenes = {}

    for c in ortho_matrix.columns:
        orfs = ortho_matrix[c]
        genes = na_model_genes
        orfs2 = orfs[orfs == 1].index.tolist()
        genes2 = genes[genes >= 50].index.tolist()
        unannotated = set(genes2) - set(orfs2)

        data = join(nucls_dir, "reactome_vs_%s.txt" % c)
        cols = [
            "gene",
            "subject",
            "PID",
            "alnLength",
            "mismatchCount",
            "gapOpenCount",
            "queryStart",
            "queryEnd",
            "subjectStart",
            "subjectEnd",
            "eVal",
            "bitScore",
        ]
        data = pd.read_csv(data, sep="\t", names=cols)
        pseudogenes[c] = {}
        unannotated_data = data[data["gene"].isin(list(unannotated))]
        for i in unannotated_data.index:
            gene = data.loc[i, "gene"]
            contig = data.loc[i, "subject"]
            start = data.loc[i, "subjectStart"]
            end = data.loc[i, "subjectEnd"]
            seq = extract_seq(join(target_genome_dir, "%s.fna" % c), contig, start, end)
            if "*" in seq:
                print(seq)
                pseudogenes[c][gene] = seq
                unannotated.discard(gene)

        print(c, unannotated)
        ortho_matrix_w_unannotated.loc[list(unannotated), c] = 1
        for g in unannotated:
            geneIDs_matrix_w_unannotated.loc[g, c] = "%s_ortholog" % g

    ortho_matrix_w_unannotated.to_csv(join(present_absence_dir, "ortho_matrix.csv"))
    geneIDs_matrix_w_unannotated.to_csv(join(present_absence_dir, "geneIDs_matrix.csv"))

    hom_matrix = pd.read_csv(join(present_absence_dir, "ortho_matrix.csv"))
    hom_matrix = hom_matrix.set_index("Unnamed: 0")

    missing_counts = {}
    initial_counts = {}
    for strain in hom_matrix.columns:
        currentStrain = hom_matrix[strain]
        nonHomologous = currentStrain[currentStrain == 0.0]
        nonHomologous = nonHomologous.index.tolist()
        for artificial in ["spontaneous", "EXCHANGE", "BIOMASS", "Diffusion", "GAP", "DEMAND", "SINK", "ORPHAN"]:
            if artificial in nonHomologous:
                nonHomologous.remove(artificial)
        missing_counts[strain] = len(nonHomologous)
        toDelete = []
        for gene in nonHomologous:
            toDelete.append(model.genes.get_by_id(gene))

        modelCopy = model.copy()
        remove_genes(modelCopy, toDelete, remove_reactions=True)
        modelCopy.id = str(strain)
        initial_model_path = join(initial_models_dir, strain + ".json")
        cobra.io.save_json_model(modelCopy, str(initial_model_path), pretty=False)
        initial_model_loaded = cobra.io.load_json_model(initial_model_path)
        initial_counts[strain] = (len(initial_model_loaded.genes), len(initial_model_loaded.reactions))

    final_counts = {}
    models = glob("%s/*.json" % initial_models_dir)
    geneIDs_matrix = pd.read_csv(join(present_absence_dir, "geneIDs_matrix.csv"))
    geneIDs_matrix = geneIDs_matrix.set_index("Unnamed: 0")

    from cobra.manipulation.modify import rename_genes
    missing_mappings = []
    for mod in models:
        model = cobra.io.load_json_model(mod)
        for column in geneIDs_matrix.columns:
            if column in mod:
                currentStrain = column
    
        IDMapping = geneIDs_matrix[currentStrain].to_dict()
        for k, v in IDMapping.items():
            if not (pd.notna(v) and v != "None" and isinstance(v, str)):
                presence = None
                if k in hom_matrix.index and currentStrain in hom_matrix.columns:
                    presence = hom_matrix.loc[k, currentStrain]
                missing_mappings.append(
                    {"strain": currentStrain, "model_gene": k, "mapped_gene": v, "present": presence}
                )
        IDMappingParsed = {k: v for k, v in IDMapping.items() if pd.notna(v) and v != "None" and isinstance(v, str)}
        rename_genes(model, IDMappingParsed)
        cobra.io.save_json_model(model, mod, pretty=False)
        # Report "initial draft" using the post-rename initial model, which matches historical behavior.
        initial_counts[currentStrain] = (len(model.genes), len(model.reactions))
    if missing_mappings:
        pd.DataFrame(missing_mappings).to_csv(join(present_absence_dir, "missing_gene_mappings.csv"), index=False)

    models = glob("%s/*.json" % initial_models_dir)
    for i in range(len(models)):
        model = cobra.io.load_json_model(models[i])
        strain = models[i].replace(initial_models_dir + "/", "")
        ort = []
        for ge in model.genes:
            if "ortholog" in str(ge):
                ort.append(ge)
        modelCopy = model.copy()
        remove_genes(modelCopy, ort, remove_reactions=True)
        modelCopy.id = str(strain)
        cobra.io.json.save_json_model(modelCopy, str(join(output_models_dir, strain + ".json")), pretty=False)
        final_counts[strain.replace(".json", "")] = (len(modelCopy.genes), len(modelCopy.reactions))

    print(f"\nPangenome reference -> {ref_gene_count} genes, {ref_rxn_count} reactions")
    print("------------------------------------------------")
    for strain in hom_matrix.columns:
        init_genes, init_rxns = initial_counts.get(strain, (None, None))
        final_genes, final_rxns = final_counts.get(strain, (None, None))
        print(f"{strain}")
        print(f"Initial draft (after removing absent genes) -> {init_genes} genes, {init_rxns} reactions")
        print(f"Final draft (after removing ortholog placeholders) -> {final_genes} genes, {final_rxns} reactions")
        print("")


if __name__ == "__main__":
    main()
    
