from step3_pdf_report import create_report

def create_pdf_inputs(output_dict, output_file):

    ################second_page//METHOD############################
        
    second_page_header = output_dict["second_page_header"] 
    second_page_text1 = output_dict["second_page_text1"]
    second_page_text2 = output_dict["second_page_text2"]

    ################Quality Check page############################

    quality_page_header = output_dict["quality_page_header"] 
    quality_page_text1 = output_dict["quality_page_text1"]
    quality_statistics = output_dict["quality_page"]
    quality_per_seq = output_dict["quality_per_seq_image_path"]


    ################Alpha diversity page############################

    alpha_diversity_page_header1 = output_dict["alpha_diversity_page_header1"]
    alpha_diversity_page_header2 = output_dict["alpha_diversity_page_header2"] 
    alpha_diversity_text1 = output_dict["alpha_diversity_page_text1"]
    alpha_diversity_table = output_dict["alpha_diversity_table"]
    alpha_diversity_text2 = output_dict["alpha_diversity_page_text2"]


    ################third_page//BIOINFORMATIC ANALYSIS#############

    third_page_header = output_dict["third_page_header"]
    third_page_text = output_dict["third_page_text"]
    third_page_phylum = output_dict["phylum_image_path"]
    third_page_class  = output_dict["class_image_path"]

    ################forth_page//images######################

    forth_page_order = output_dict["order_image_path"]
    forth_page_family = output_dict["family_image_path"]

    #############fifth_page///images##############################
    fifth_page_genus = output_dict["genus_image_path"]
    fifth_page_species = output_dict["species_image_path"]

    ###################sixth_page//raw_data########################

    raw_data = output_dict["all_summary"]

    #####################sample_name##############################

    sample_name = output_dict["sample_name"]


    report_content = {
                "second_page": {"header":second_page_header, "text1": second_page_text1, "text2":second_page_text2},
                "quality_page": {"header":quality_page_header, "text1": quality_page_text1, "quality_statistics":quality_statistics, "quality_per_seq_image_path":quality_per_seq},
                "alpha_diversity_page": {"header1":alpha_diversity_page_header1, "header2":alpha_diversity_page_header2, "text1": alpha_diversity_text1, "alpha_diversity_table": alpha_diversity_table, "text2":alpha_diversity_text2},
                "third_page": {"header": third_page_header, "text": third_page_text, "phylum":third_page_phylum, "class": third_page_class},
                "forth_page": {"order": forth_page_order, "family": forth_page_family}, 
                "fifth_page": {"genus": fifth_page_genus, "species": fifth_page_species}, 
                "raw_data": raw_data, "sample_name":sample_name
    }

    create_report(report_content, output_file)