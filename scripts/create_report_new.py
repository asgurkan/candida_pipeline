import sys

import pandas as pd
import my_plots as mp

# inputs
INPUT_SUMMARY = snakemake.input.summary
INPUT_ALPHA = snakemake.input.alpha
INPUT_QUALITY_DATA = snakemake.input.quality_per_seq

# outputs
OUTPUT_EXCEL = snakemake.output.excel 

# params
SAMPLE_NAME = snakemake.params.sample_name
CSS_STYLE = snakemake.params.css_style
ORGANISM = snakemake.params.organism

# functions
def detect_domain_count(df, sample_name):
    for i, data in df.iterrows():
        if data["taxon"] == ORGANISM:
            domain_count = df.loc[i,str(sample_name)]
            break
    return domain_count 


project_name = OUTPUT_EXCEL.split("/")[-2]

OUTPUT_FILES_DICT = dict()

# Excel oluşturma; input csv formatına dönüştürülüp, taxon columnundaki "root"lar çıkarılıyor.
summary = pd.read_csv(INPUT_SUMMARY, sep="\t")
alpha_summary = pd.read_csv(INPUT_ALPHA, sep="\t")

# Root olan satırları çıkar
summary = summary[summary["taxon"] != "Root"]

# tax_level dictionary'si oluşturularak her bir sayının uygun tax level'ı giriliyor.
tax_levels = {1: "Domain", 2: "Phylum", 3: "Class", 4: "Order", 5: "Family", 6: "Genus", 7: "Species"}
summary["Taxonomy"] = summary["taxlevel"].map(lambda x: tax_levels[x])

# İşe yaramayacak bazı sütunlar siliniyor ve index resetleniyor. rankID column'un ismi MassID olarak değiştiriliyor.
summary.drop(columns=["daughterlevels", "taxlevel", "total"], inplace=True)
summary.reset_index(inplace=True)
summary.rename(columns={"rankID": "MassID"}, inplace=True)

# Summary DataFrame'inin kopyası alınarak domain_count hesaplanıyor.
summary_1 = summary.copy()
print(summary_1)
domain_count = detect_domain_count(summary_1, SAMPLE_NAME)

# DataFrame'e bir abundance sütunu eklenerek yüzde üzerinden hesaplamalar yapılıyor.
summary_1["abundance"] = summary_1[SAMPLE_NAME].apply(lambda x: ((int(x) / int(domain_count)) * 100))

# Örnek ismi sütunu değiştirilip, abundance sütunu olarak yeniden adlandırılıyor ve sıralanıyor.
summary_1.rename(columns={SAMPLE_NAME: f"{SAMPLE_NAME}_Count"}, inplace=True)
summary_1.rename(columns={"abundance": SAMPLE_NAME}, inplace=True)
summary_1.sort_values(by=SAMPLE_NAME, ascending=False, inplace=True)

# Sütun isimleri ve sıralamaları düzenleniyor.
s_summary = summary_1.copy()
s_summary.rename(columns={SAMPLE_NAME: f"%{SAMPLE_NAME}"}, inplace=True)
s_summary.rename(columns={f"{SAMPLE_NAME}_Count": SAMPLE_NAME}, inplace=True)
s_summary = s_summary[["index", "MassID", "taxon", "Taxonomy", SAMPLE_NAME, f"%{SAMPLE_NAME}"]]

# Farklı sheetler için tablolar oluşturulacak.
s_summary_sheet1 = s_summary.copy()

try:
    with pd.ExcelWriter(OUTPUT_EXCEL) as writer:
        # Summary sheetini excele yazılıyor.
        s_summary_sheet1.sort_values(by=SAMPLE_NAME, ascending=False, inplace=True)
        OUTPUT_FILES_DICT["all_summary"] = s_summary_sheet1
        s_summary_sheet1.to_excel(writer, sheet_name="Summary", index=False)

        # Species sheetini excele yazılıyor.
        s_summary_species = s_summary_sheet1.query('Taxonomy == "Species"')
        s_summary_species = s_summary_species.drop(["index", "MassID", "Taxonomy"], axis=1)
        s_summary_species.to_excel(writer, sheet_name='Species', index=False)
        OUTPUT_FILES_DICT["species"] = s_summary_species

        # Genus sheetini excele yazılıyor.
        s_summary_genus = s_summary_sheet1.query('Taxonomy == "Genus"')
        s_summary_genus = s_summary_genus.drop(["index", "MassID", "Taxonomy"], axis=1)
        s_summary_genus.to_excel(writer, sheet_name='Genus', index=False)
        OUTPUT_FILES_DICT["genus"] = s_summary_genus

        # Family sheetini excele yazılıyor.
        s_summary_family = s_summary_sheet1.query('Taxonomy == "Family"')
        s_summary_family = s_summary_family.drop(["index", "MassID", "Taxonomy"], axis=1)
        s_summary_family.to_excel(writer, sheet_name='Family', index=False)
        OUTPUT_FILES_DICT["family"] = s_summary_family

        # Order sheetini excele yazılıyor.
        s_summary_order = s_summary_sheet1.query('Taxonomy == "Order"')
        s_summary_order = s_summary_order.drop(["index", "MassID", "Taxonomy"], axis=1)
        s_summary_order.to_excel(writer, sheet_name="Order", index=False)
        OUTPUT_FILES_DICT["order"] = s_summary_order

        # Phylum sheetini excele yazılıyor.
        s_summary_phylum = s_summary_sheet1.query('Taxonomy == "Phylum"')
        s_summary_phylum = s_summary_phylum.drop(["index", "MassID", "Taxonomy"], axis=1)
        s_summary_phylum.to_excel(writer, sheet_name='Phylum', index=True)
        OUTPUT_FILES_DICT["pyhlum"] = s_summary_phylum


        # Alpha diversity sheetini excele yazılıyor.
        alpha_summary = alpha_summary.drop(["chao1_ci"],axis=1)
        alpha_column_names = ["Sample", "Shannon", "Chao1", "Simpson", "Enspie", "Fisher Alpha", "Mcintosh D", "Mcintosh E"] 
        alpha_summary.columns = alpha_column_names
        alpha_summary.to_excel(writer, sheet_name="Alpha_Diversity", index=False)

        # Pie chart oluşturuluyor.
        level_dict_piechart = mp.create_pie_charts(s_summary, OUTPUT_EXCEL, SAMPLE_NAME)
        values_level_dict = level_dict_piechart.values()

        for a in values_level_dict:
            worksheet = writer.sheets['Summary']
            worksheet.insert_image("I2", a)

        # Per seq quality chart oluşturuluyor
        input_quality_data_df = pd.read_csv(INPUT_QUALITY_DATA, sep="\t")
        quality_png, quality_html = mp.create_per_seq_quality_chart(input_quality_data_df, OUTPUT_EXCEL, SAMPLE_NAME)

except Exception as e:
    print(f"An error occurred: {e}")
    sys.exit(1)
