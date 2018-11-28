#!/usr/bin/env python
#title           :functions.py
#description     :Functions and methods for SMAC
#author          :s.pirro
#date            :20180123
#version         :0.1
#notes           :
#==============================================================================

# importing libraries
import os, sys, glob, subprocess, re, logging, urllib
from Bio import Entrez, Medline
from collections import Counter
from dateutil import parser # library to convert any type of text into datetime format
from vars import * # importing custom variables

## functions to split array in fixed chunks
def split(arr, size):
     arrs = []
     while len(arr) > size:
         pice = arr[:size]
         arrs.append(pice)
         arr   = arr[size:]
     arrs.append(arr)
     return arrs

## Function to retrieve papers from PubMed
def RetrievePapersFromPubmed(term, limit, email):
    Entrez.email = email
    pubmed_ids = []
    # retrieving the total number of papers
    total_pubs = Entrez.esearch(db='pubmed', term=term, sort="relevance")
    result = Entrez.read(total_pubs)
    max_count = list(range(0, int(result["Count"])))
    print('Total number of publications containing {0}: {1}'.format(term, result['Count']))

    # we will perform multiple queries by chunks (to speedup process)
    chunks = split(max_count, 10000)
    # check if limit is a number or not
    if limit == "all":
        for chunk in chunks:
            total_pubs = Entrez.esearch(db='pubmed', retmax=len(chunk), retstart=chunk[0], term=term)
            # reading results
            result = Entrez.read(total_pubs)
            pubmed_ids.extend(result['IdList'])
    else:
        limit_array = list(range(0,int(limit)))
        chunks = split(limit_array, 10000)
        for chunk in chunks:
            total_pubs = Entrez.esearch(db='pubmed', retmax=len(chunk), retstart=chunk[0], term=term)
            # reading results
            result = Entrez.read(total_pubs)
            pubmed_ids.extend(result['IdList'])

    return pubmed_ids

## Fuction to retrieve infos from pmids
def RetrieveInfoFromPMID(pubmed_ids, data_folder, email):
    ### initialising email for the queries
    Entrez.email = email

    ## retrieve GSE (GEO series ids), platforms and links from selected publication
    print("## Retrieving GEO datasets for the papers...(this may take long)")
    all_geo_info = RetrieveGeoFromPMID(pubmed_ids)

    print("\n## Retrieving info for the papers...")

    # here we use the epost util to speedup process
    search_handle = Entrez.epost(db='pubmed', id=",".join(pubmed_ids))
    search_results = Entrez.read(search_handle)
    # Getting web environment and query_key (for efetch after)
    webenv, query_key = search_results["WebEnv"], search_results["QueryKey"]
    # Initialising papers info dictionary and the mesh heading container (list)
    papers_infos = {}
    mesh_headings_container = []
    # This analysis will be performed by chuncks of 10000 pmids each (for not overloading the servers)
    chunks = split(pubmed_ids, 10000)
    # Processing the epost and efetch process for each chunk
    for i,chunk in enumerate(chunks):
        PrintInline("Parsing {0} chunk (at least 10,000) of papers...".format(i+1))
        # we can perform multiple efetch calls, remembering the index (index * 1,000)
        pmids_info = Entrez.efetch(db='pubmed', rettype='medline', retmode='text', retstart=i*1000, retmax=len(chunk), webenv=webenv, query_key=query_key)
        records = Medline.parse(pmids_info)

        ## iterating into records set
        for record in records:
            # Getting vars and assign defualt if not exist
            pmid = record.get('PMID', '')
            au = record.get('AU', '') #author
            journal = record.get('JT', '') # journal title
            title = record.get('TI', '') #title
            abstract = record.get('AB', '') # abstract
            # in case of strage dates, we add 1st January 2018 as default
            try:
                pub_date = parser.parse(record.get('DP', '')).strftime("%Y-%m-%d")
            except ValueError:
                pub_date = parser.parse("2018-01-01").strftime("%Y-%m-%d")
            mesh_headings = record.get('MH', '') # mesh terms

            if pmid in all_geo_info.keys():
                gse_codes = all_geo_info[pmid]["gse_codes"]
                platforms = all_geo_info[pmid]["platforms"]
                ftp_links = all_geo_info[pmid]["ftp_links"]
            else:
                gse_codes = platforms = ftp_links = ""

            # updating dictionary
            papers_infos[pmid] = {"author":au, "journal":journal, "title":title, "abstract":abstract, "pub_date":pub_date, "mesh_headings":mesh_headings, "gse_codes":gse_codes, "platforms":platforms, "ftp_links":ftp_links, "analysis":"gene_expression,correlation,gene_networks"}
            # updating mesh heading list
            mesh_headings_container.extend(mesh_headings)

    # replacing * in mesh headings
    mesh_headings_container_no_unique = [mh.replace('*', '') for mh in mesh_headings_container]
    mesh_headings_container = list(set(mesh_headings_container_no_unique))
    return pubmed_ids, mesh_headings_container_no_unique, mesh_headings_container, papers_infos

## function to retrieve geo info from pmid
def RetrieveGeoFromPMID(pmid):
    # we have to split the pmids for limitations by NCBI
    chunks = split(pmid, 1000)
    geo_info_container = {}
    for i,chunk in enumerate(chunks):
        PrintInline("Parsing {0} chunk (at least 1,000) of papers...".format(i+1))
        try:
            geo_query = Entrez.elink(dbfrom='pubmed', db='gds', id=list(map(int, chunk)))
            geo_info = Entrez.read(geo_query)
        # if the query fails for some reasons, try to reconnect to Pubmed after 5 seconds
        except (urllib.HTTPError, urllib.URLError, RuntimeError, KeyError) as e:
            WaitForSucceed(e)
        for gi in geo_info:
            gds_codes = []
            pmid = gi["IdList"][0]
            if (len(gi["LinkSetDb"]) > 0): # this means we found some results
                for gds in gi["LinkSetDb"][0]["Link"]:
                    gds_codes.append(gds["Id"])

                # retrieve gse codes and other informations
                gds_codes = ",".join(list(set(gds_codes)))
                geo_query = Entrez.esummary(db='gds', id=gds_codes)
                geo_info = Entrez.read(geo_query)
                geo_query.close()
                # initialising arrays
                gse_codes, platforms, ftp_links = ([] for i in range(3))
                # iterating retrieved output and extract infos
                for i in range(0, len(geo_info)):
                    gse_codes.append(geo_info[i].get("GSE", "NA"))
                    platforms.append(geo_info[i].get("GPL", "NA"))
                    ftp_links.append(geo_info[i].get("FTPLink", "NA"))

                # uniquing arrays and returning
                gse_codes = list(set(gse_codes))
                platforms = list(set(platforms))
                ftp_links = list(set(ftp_links))
                geo_info_container[pmid] = {"gse_codes": gse_codes, "platforms":platforms, "ftp_links":ftp_links}
            else:
                next

    return geo_info_container

## function to retrieve gene list from articles
def RetrieveGenesFromPMID(pmid, data_folder):
    geo_query = Entrez.elink(dbfrom='pubmed', db='gene', cmd='neighbor_score', id=pmid)
    genes_info_container = {}
    try:
        geo_info = Entrez.read(geo_query)
        for gi in geo_info:
            try:
                for gene in gi["LinkSetDb"][0]["Link"]:
                    if int(gene["Score"]) > 0:
                        genes_info_container[gene["Id"]] = gene["Score"]
            except: # there are no genes to convert
                next
    except:
        next

    # mapping genes to uniprot -- just if there are genes to convert!
    if len(genes_info_container) > 0:
        try:
            geo_query = Entrez.esummary(db='gene', id=",".join(genes_info_container.keys()))
            geo_info = Entrez.read(geo_query)
            geo_query.close()
            converted_genes = {}
            for gi in geo_info["DocumentSummarySet"]["DocumentSummary"]:
                converted_genes[gi.attributes["uid"]] = gi["Name"]

            # writing list into a file
            file_io = open("{0}/gene_list.txt".format(data_folder), "w")
            for gene, score in genes_info_container.items():
                if gene in converted_genes.keys():
                    file_io.write("{0}\t{1}\n".format(converted_genes[gene], score))
            file_io.close
        # in case we have problems in the connection...
        except (urllib.HTTPError, urllib.URLError, RuntimeError, KeyError) as e:
            WaitForSucceed()

    return genes_info_container

## function for sorting MeSH headings
def SortMesh(type, mesh_terms, data_folder):
    if type == "level":
        # open mesh dictionary
        file_io = open("src/mesh_dict.txt", "r").read().split("\n")
        # filling mesh dictionary
        mesh_dict = {}
        for line in file_io:
            try:
                parts = line.split("\t")
                if parts[1] is not "":
                    mesh_dict[parts[0]] = parts[1] # parts[0] = mh, parts[1] = level
            except IndexError:
                next

        # mapping mesh terms to the dictionary
        mesh_list = {}
        for mh in mesh_terms:
            if mh in mesh_dict.keys():
                mesh_list[mh] = int(mesh_dict[mh])

        # save mesh_list
        io_stream = open("{0}/mesh_sorted_by_level.txt".format(data_folder), "w")
        for mh, level in mesh_list.items():
            io_stream.write("{0}\t{1}\n".format(mh, level))
        io_stream.close()

        return "{0}/mesh_sorted_by_level.txt".format(data_folder)
    elif type == "laboundance":
        # count number occurrences for mesh terms
        mesh_list = Counter(mesh_terms)

        # save mesh_list
        io_stream = open("{0}/mesh_sorted_by_laboundance.txt".format(data_folder), "w")
        for mh, level in mesh_list.items():
            io_stream.write("{0}\t{1}\n".format(mh, level))
        io_stream.close()

        return "{0}/mesh_sorted_by_laboundance.txt".format(data_folder)
    elif type == "gaboundance":
        mesh_list = {}
        # count number of total papers related to the MeSH heading
        for index, mh in enumerate(mesh_terms):
            PrintInline("Retrieving global aboundance for: {0}.. -- {1}/{2}".format(mh[0:10], index+1, len(mesh_terms)))
            try:
                query = Entrez.esearch(db='pubmed', term="{0} [MeSH Terms]".format(mh), rettype='Count')
                result = Entrez.read(query)
                total_pubs = result['Count']
                mesh_list[mh] = total_pubs
            except (urllib.HTTPError, urllib.URLError, RuntimeError, KeyError) as e:
                WaitForSucceed(e)

        # save mesh_list
        io_stream = open("{0}/mesh_sorted_by_gaboundance.txt".format(data_folder), "w")
        for mh, level in mesh_list.items():
            io_stream.write("{0}\t{1}\n".format(mh, level))
        io_stream.close()

        return "{0}/mesh_sorted_by_gaboundance.txt".format(data_folder)

# this function:
#   - downloads the GEO Dataset (Using GEOQuery R),
#   - map the probes
#   - export the converted and reduced dataset in a text file
#   -
def DownloadGEODataset(gse, gse_folder):
    # Getting log of analysed files
    try:
        subprocess.check_call(["Rscript","scripts/GEO/DownloadGSE.R","--gse",gse,"--dir",gse_folder], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    # manage exec error
    except subprocess.CalledProcessError:
        print(exec_3_err % gse)
        os.rmdir(gse_folder)

# This function creates the target file from the pData downloaded from GEO
def CreateTargetFile(gse_folder):
    # listing target file
    pData_files = glob.glob("{0}/pData.*.tsv".format(gse_folder))
    # iterating into the pData files and retrieve the header (for SAMPLES and description)
    for i in range(0,len(pData_files)):
        pData_file = "{0}/pData.{1}.tsv".format(gse_folder,str(i+1))
        try:
            subprocess.check_call(["Rscript","scripts/DataAnalysis/create_target.R","-p",pData_file,"-d",gse_folder,"-c",str(i+1)], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        # manage exec error
        except subprocess.CalledProcessError:
            print(exec_5_err % "target creation")

# this function perform multiple analyses on the gse dataset, according to what chosen by the user
def PerformAnalyses(gse_folder, tFile, eFile, analyses, cont):
    # defining dictionary of commands for the selected analyses
    analysis_command = {
        "pca": {"script_name":"PCA.R", "exp": eFile, "target": tFile, "output": gse_folder},
        "receptor_status": {"script_name":"mclust.R", "exp": eFile, "target": tFile, "output": gse_folder},
        "tumour_purity": {"script_name":"estimate.R", "exp": eFile, "target": tFile, "output": gse_folder},
        "molecular_classification": {"script_name":"pam50.R", "exp": eFile, "target": tFile, "output": gse_folder}
	"gene_expression": {"script_name":"gene_expr.R", "exp": eFile, "target": tFile, "output": gse_folder}
    }

    # splitting analyses string into a list
    analyses = analyses.rstrip().split(",")

    # iterating into the analyses list and launch commands according to the dictionary
    for an in analyses:
        # check if we digit a valid analysis
        if an in analysis_command.keys():
            # collecting variables for launching the analysis
            script_name = analysis_command[an]["script_name"]
            exp_file = analysis_command[an]["exp"]
            target_file = analysis_command[an]["target"]
            output_dir = analysis_command[an]["output"]
            # Getting log of analysed files
            try:
                subprocess.check_call(["Rscript","scripts/DataAnalysis/{0}".format(script_name),"-e",exp_file,"-t",target_file,"-c","target","-d",output_dir,"-n",cont], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                PrintInline("SUCCESS -- {0}".format(an))
            # manage exec error
            except subprocess.CalledProcessError:
                print(exec_5_err % an)
        else:
            print(exec_4_err % an)

# This function prints messages inLine
def PrintInline(text):
    print("    {0}".format(text), end="\r", flush=True)

# This function update Infos for downloaded papers with performed analyses
def UpdatePapersInfo(papers_info, pmid, gse, gse_folder, analyses, data_folder):
    # initialising array with performed analyses
    perf_analyses = []
    # splitting analyses string into a list
    analyses = analyses.rstrip().split(",")
    # iterating into the analyses list and search for corrensponding files
    for an in analyses:
        n_plots = len(glob.glob("{0}/{1}*.html".format(gse_folder, an)))
        perf_analyses.append(an) if n_plots > 0 else False

    # update file with performed analyses (by GSE and PMID)
    file_io = open("{0}/performed_analyses.tsv".format(data_folder), "a")
    file_io.write("{0}\tGSE{1}\tgene_expression,correlation,gene_networks,{2}\n".format(pmid, gse, ",".join(perf_analyses)))
    file_io.close()
    # update papers_info with the performed analysis
    papers_info[pmid]["analysis"] += ",{0}".format(",".join(perf_analyses))

    return papers_info

# This function creates a final literature report in a file
def CreateLiteratureReport(papers_info, data_folder):
    # creating file
    file_io = open("{0}/literature.tsv".format(data_folder), "w")
    # writing header to the file
    file_io.write("PMID\tTitle\tJournal\tPubDate\tAuthors\tAbstract\tMeSH\tAnalysis\tCurated\n")
    for pmid, infos in papers_info.items():
        # setting default values
        title = infos.get("title", "NA")
        journal = infos.get("journal", "NA")
        pub_date = infos.get("pub_date", "NA")
        authors = infos.get("author", "NA")
        abstract = infos.get("abstract", "NA")
        mesh = infos.get("mesh_headings", "NA")
        analysis = infos.get("analysis", "NA")
        file_io.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(pmid,title,journal,pub_date,",".join(authors),abstract,",".join(mesh),analysis,"No"))
    # close file
    file_io.close()

# This function cleans samples definitions from "background terms"
def CleanDef(background_list, definition):
    # copying the text in a new var (with all lower cases)
    new_definition = definition.lower()
    # iterate into regex in the list and sub what discovered with empty chars
    for rgx in background_list:
        new_definition = re.sub(rgx, '', new_definition)
    # returning the new definitions
    return new_definition

# This function wait 5 seconds if the connection gives any error
def WaitForSucceed(error):
    # take not of the error
    logging.error(error)
    # debugging the error by waiting for 5 second
    logging.debug("wait 5 seconds to reconnect...")
    time.sleep(5)
