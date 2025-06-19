import pandas as pd
from step2_create_pdf_file import create_pdf_inputs

# inputs
SUMMARY = snakemake.input.summary
BEFORE_PER_SEQ_QUALITY_SCORES = snakemake.input.before_per_seq_quality_scores
BEFORE_GENERAL_STATS = snakemake.input.before_general_stats
AFTER_PER_SEQ_QUALITY_SCORES = snakemake.input.after_per_seq_quality_scores
AFTER_GENERAL_STATS = snakemake.input.after_general_stats
ALPHA_DIVERSITY = snakemake.input.alpha


# outputs
OUTPUT_PDF = snakemake.output.pdf 

# params
SAMPLE_NAME = snakemake.params.sample_name

output_files_dict = dict()

#input csv formatına dönüştürülüp, taxon columnundaki "root"lar çıkarılıyor.  
summary = pd.read_csv(SUMMARY, sep="\t")

summary = summary[summary["taxon"]!="Root"]
tax_levels = {1:"Domain", 2:"Phylum", 3:"Class", 4:"Order", 5:"Family", 6:"Genus", 7:"Species"}

summary["Taxonomy"] = summary["taxlevel"].map(lambda x: tax_levels[x])
summary.drop(columns=["daughterlevels","taxlevel","total"], inplace=True)
summary.reset_index(inplace=True)
summary.rename(columns={"rankID":"MassID"}, inplace=True)

output_files_dict["sample_name"] = SAMPLE_NAME 

# page six/ summary table
summary_1 = summary.copy()
summary_1 = summary_1[summary_1["Taxonomy"] == "Species"]
species_count_df = summary_1[summary_1["Taxonomy"] == "Species"] #### DOMAIN'den SPECIES'E ÇEVİRDİM. -EA 
species_count = species_count_df[SAMPLE_NAME].sum()
summary_1["Abundance (%)"] = summary_1[SAMPLE_NAME].apply(lambda x: (x/species_count*100).round(3))
summary_1.rename(columns = {"taxon": "Taxon"}, inplace=True)
#ummary_1.rename(columns = {"taxon": "Takson"}, inplace=True)

summary_1.rename(columns={SAMPLE_NAME: "Read Count"}, inplace=True)
#summary_1.rename(columns={SAMPLE_NAME: "Okuma Sayısı"}, inplace=True)
#summary_1.rename(columns = {"Taxonomy": "Taksonomi"}, inplace=True)
#summary_1.rename(columns = {"Abundance (%)": "Bolluk (%)"}, inplace=True)

summary_1 = summary_1[["Taxon","Taxonomy", "Read Count", "Abundance (%)"]] 
#summary_1 = summary_1[["Takson","Taksonomi", "Okuma Sayısı", "Bolluk (%)"]] 
summary_1 = summary_1.sort_values("Abundance (%)", ascending=False)
output_files_dict["all_summary"] = summary_1

# Quality table
before_general_stats = pd.read_csv(BEFORE_GENERAL_STATS, sep="\t")
after_general_stats = pd.read_csv(AFTER_GENERAL_STATS, sep="\t")

before_general_stats.columns = [f"{i}_before_trim" for i in before_general_stats.columns]
after_general_stats.columns = [f"{i}_after_trim" for i in after_general_stats.columns]

after_general_stats = after_general_stats.drop(["Measure_after_trim"], axis=1)
general_stats = pd.concat([before_general_stats, after_general_stats], axis=1)
general_stats.columns = ["Quality Specifications", "Quality Values Before Trim", "Quality Values After Trim"]
#general_stats.columns = ["Kalite Özellikleri", "Trimden Önce Kalite Değerleri", "Trimden Sonra Kalite Değerleri"]

output_files_dict["quality_page"] = general_stats

# Alpha diversity
alpha_diversity = pd.read_csv(ALPHA_DIVERSITY, sep="\t")
alpha_diversity = alpha_diversity.drop(["chao1_ci"],axis=1)
alpha_column_names = ["Sample", "Shannon", "Chao1", "Simpson", "Enspie", "Fisher Alpha", "Mcintosh D", "Mcintosh E"] 
alpha_diversity.columns = alpha_column_names
alpha_column_names_round = ["Shannon", "Simpson", "Enspie", "Fisher Alpha", "Mcintosh D", "Mcintosh E"] 
alpha_diversity[alpha_column_names_round] = alpha_diversity[alpha_column_names_round].apply(lambda x: round(x,3))
output_files_dict["alpha_diversity_table"] = alpha_diversity

# charts
phylum_image_path = snakemake.params.phylum
class_image_path = snakemake.params.class_pic
order_image_path = snakemake.params.order
family_image_path = snakemake.params.family
genus_image_path = snakemake.params.genus
species_image_path = snakemake.params.species
quality_per_seq_image_path   = snakemake.params.quality


output_files_dict["phylum_image_path"] = phylum_image_path
output_files_dict["class_image_path"] = class_image_path
output_files_dict["order_image_path"] = order_image_path
output_files_dict["family_image_path"] = family_image_path
output_files_dict["genus_image_path"] = genus_image_path
output_files_dict["species_image_path"] = species_image_path
output_files_dict["quality_per_seq_image_path"] = quality_per_seq_image_path

output_files_dict["second_page_header"] = "1. METHOD"
#output_files_dict["second_page_header"] = "1. METOD"
output_files_dict["second_page_text1"] = ("For the library preparation of the sample, Ligation sequencing kit (SQK-NBD114.24;" 
                                        "Oxford Nanopore Technologies) protocol was used. First the sample was end-prepped,"
                                        "dA tailed and nick-repaired in a final reaction volume of 60 μl, then the sample was purified using MobiomX MagBeads "
                                        "and adapters were ligated to the repaired ends. Samples were then purified "
                                        "and quantified spectroflorimetrically. Total library was loaded onto a Spot-On flow cell (FLO-MIN114) after priming "
                                        "and the sequencing run was initiated on an Mk1C™ device (Oxford Nanopore Technologies) using the latest MinKNOW™ software. "
                                        "Sequencing was stopped when enough data was obtained or the maximum run time of 72 hours was completed.")
#output_files_dict["second_page_text1"] = ("Numunenin kütüphane hazırlığı için Ligation sequencing kit (SQK-NBD114.24;” Oxford Nanopore Technologies) protokolü kullanılmıştır. İlk olarak numune 60 μl'lik nihai reaksiyon hacminde uç hazırlandı, dA kuyruklandı ve nik onarıldı, ardından numune MobiomX MagBeads kullanılarak saflaştırıldı ve adaptörler onarılan uçlara bağlandı. Örnekler daha sonra saflaştırıldı ve spektroflorimetrik olarak ölçüldü. Toplam kütüphane, priming işleminden sonra bir Spot-On akış hücresine (FLO-MIN114) yüklendi ve sekanslama çalışması en son MinKNOW™ yazılımı kullanılarak bir Mk1C™ cihazında (Oxford Nanopore Technologies) başlatıldı. Yeterli veri elde edildiğinde veya 72 saatlik maksimum çalışma süresi tamamlandığında sekanslama durduruldu.")

output_files_dict["second_page_text2"] = ("After sequencing, the results obtained in pod5 format were converted to fastq format using latest version of dorado software "
                                          "(base-calling and de-multiplexing). Barcode and adapter sequences were cleared using dorado software, and universal primers and "
                                          "tags were also deleted by deleting 15 bases from both ends of the sequences. After clearing the sequences, reads 1250-1750 bp "
                                          "long were filtered with Trimmomatic, and the remaining reads were excluded from the analysis. "
                                          "The trimmed reads were aligned to the Massive Bioinformatics database produced using 16S rRNA genes of bacterial whole genomes submitted to "
                                          "NCBI using the Massive Align Engine (MAE). Using the taxonomic details of the 16S rRNA genes to which the reads were aligned and the number of "
                                          "reads aligned to the genes, a biome file was produced. Using the capabilities offered by the qiime2 platform, alpha diversity analysis was "
                                          "carried out using several kinds of indexes in order to conduct phylogenetic analysis using the biome file which had been produced. "
                                          "Utilizing the Mothur platform, taxonomic classifications were arranged. The Python programming language's libraries were used to create the "
                                          "tables and graphs in the analysis.")

# output_files_dict["second_page_text2"] = ("Sekanslamanın devamında, pod5 formatında elde edilen sonuçlar dorado yazılımının son versiyonu kullanılarak fastq formatına dönüştürüldü (base-calling ve de-multiplexing). Barkod ve adaptör dizileri dorado yazılımı kullanılarak temizlendi ve evrensel primerler ve etiketler de dizilerin her iki ucundan 15 baz silinerek silindi. Diziler temizlendikten sonra, 1250-1750 bp uzunluğundaki okumalar Trimmomatic ile filtrelendi ve kalan okumalar analizden çıkarıldı. Filtrelenmiş okumalar, Massive Align Engine (MAE) kullanılarak NCBI’da bulunan bakteriyel tüm genomların 16S rRNA genleri elde edilerek üretilen Massive Bioinformatics’e ait 16S rRNA veritabanına hizalandı. Okumaların hizalandığı 16S rRNA genlerinin taksonomik detayları ve genlere hizalanan okuma sayısı kullanılarak bir biom dosyası üretildi. Oluşturulan biom dosyasını kullanarak filogenetik analiz yapmak için qiime2 platformunun sunduğu analizler kullanıldı. Çeşitli indeksler kullanılarak alfa çeşitlilik analizi gerçekleştirildi. Mothur platformu kullanılarak taksonomik sınıflandırmalar düzenlendi. Analizde tablo ve grafiklerin oluşturulması için Python programlama dilinin kütüphaneleri kullanıldı.")


output_files_dict["quality_page_header"] = "2. General Quality Statistics"
# output_files_dict["quality_page_header"] = "2. Genel Kalite İstatistikleri"


output_files_dict["quality_page_text1"] = ("Quality analysis was performed on the sequencing data before and after trimming using FastQC."
                                           "The tables below provide an overview of the quality metrics, including total sequences, total bases, "
                                           "sequence length, and GC content. These metrics are important for assessing the integrity and reliability "
                                           "of the sequencing data.")
# output_files_dict["quality_page_text1"] = ("Sekanslama verilerinin, trimlenme işleminden önce ve sonra FastQC kullanılarak kalite analizi yapılmıştır. Aşağıdaki tablolar, toplam sekans sayısı, toplam baz sayısı, sekans uzunlukları ve GC içeriği dahil olmak üzere kalite ölçütlerine genel bir bakış sunmaktadır. Bu metrikler, sekanslama verilerinin bütünlüğünü ve güvenilirliğini değerlendirmek için önemlidir.")


output_files_dict["alpha_diversity_page_header1"] = "3. BIOINFORMATICS ANALYSIS"
# output_files_dict["alpha_diversity_page_header1"] = "3. BIOINFORMATİK ANALİZLER"


output_files_dict["alpha_diversity_page_header2"] = "3.1 Alpha Diversity"
# output_files_dict["alpha_diversity_page_header2"] = "3.1 Alfa Çeşitlilik"

output_files_dict["alpha_diversity_page_text1"] = ("Refers to the diversity within a specific area or ecosystem, "
                                                   "often measured by the number of species (species richness) and their "
                                                   "relative abundances. It provides a snapshot of the diversity at a single "
                                                   "location or in a single sample, capturing the complexity and variety of "
                                                   "species present. Alpha diversity helps to assess the health and stability of "
                                                   "an ecosystem or community.")
# output_files_dict["alpha_diversity_page_text1"] = ("Alfa çeşitlilik ,belirli bir alan veya ekosistemdeki çeşitliliği ifade eder, genellikle tür sayısı (tür zenginliği) ve göreceli bollukları ile ölçülür. Mevcut türlerin karmaşıklığını ve çeşitliliğini yakalayarak tek bir konumdaki veya tek bir örnekteki çeşitliliğin anlık görüntüsünü sağlar. Alfa çeşitliliği, bir ekosistemin veya topluluğun sağlığını ve istikrarını değerlendirmeye yardımcı olur.")


output_files_dict["alpha_diversity_page_text2"] = (
                                                    "\n**Shannon:** Measures the diversity of a community by taking into account both the number of species and "
                                                    "their relative abundance. It is based on the concept of entropy from information theory.\n\n"
                                                    "**Chao1:** An estimator of species richness that predicts the number of species in a community based on "
                                                    "the number of rare species and the number of singletons and doubletons.\n\n"
                                                    "**Simpson:** Quantifies the diversity of a community by considering the probability that two randomly selected "
                                                    "individuals belong to the same species. It emphasizes the dominance of species.\n\n"
                                                    "**Evenness (Enspie):** A measure of how equal the abundances of different species are in a community. It reflects "
                                                    "how evenly the individuals are distributed among the species.\n\n"
                                                    "**Fisher Alpha**: An index used to estimate species richness based on the frequency of species occurrence. It is part "
                                                    "of the family of abundance-based diversity measures.\n\n"
                                                    "**McIntosh D:** A diversity index that accounts for both the number of species and their evenness, providing a measure "
                                                    "of biodiversity that balances species richness and equitability.\n\n"
                                                    "**McIntosh E:** An evenness index that complements McIntosh D by focusing specifically on the evenness of species "
                                                    "distribution, helping to understand the relative abundance of species in a community."
                                                )

# output_files_dict["alpha_diversity_page_text2"] = (
#                                                     "\n**Shannon:** Hem türlerin sayısını hem de göreceli bolluklarını dikkate alarak bir topluluğun çeşitliliğini ölçer. Enformasyon teorisindeki entropi kavramına dayanır.\n\n"
#                                                     "**Chao1:** Nadir türlerin sayısına ve tekil ve ikili türlerin sayısına dayalı olarak bir topluluktaki tür sayısını tahmin eden bir tür zenginliği tahmincisidir.\n\n"
#                                                     "**Simpson:** Rastgele seçilen iki bireyin aynı türe ait olma olasılığını göz önünde bulundurarak bir topluluğun çeşitliliğini ölçer. Türlerin baskınlığını vurgular.\n\n"
#                                                     "**Evenness (Enspie):** Bir toplulukta farklı türlerin bolluklarının ne kadar eşit olduğunun bir ölçüsüdür. Bireylerin türler arasında ne kadar eşit dağıldığını yansıtır.\n\n"
#                                                     "**Fisher Alpha**: Türlerin ortaya çıkma sıklığına dayalı olarak tür zenginliğini tahmin etmek için kullanılan bir indeks. Bolluğa dayalı çeşitlilik ölçümleri ailesinin bir parçasıdır.\n\n"
#                                                     "**McIntosh D:** Tür zenginliği ve eşitliğini dengeleyen bir biyoçeşitlilik ölçüsü sağlayan, hem tür sayısını hem de türlerin eşitliğini hesaba katan bir çeşitlilik endeksidir.\n\n"
#                                                     "**McIntosh E:** Özellikle tür dağılımının düzgünlüğüne odaklanarak McIntosh D'yi tamamlayan ve bir topluluktaki türlerin göreceli bolluğunu anlamaya yardımcı olan bir evenness endeksidir."
#                                                 )


output_files_dict["third_page_header"] = "3.2 Abundance Information According to Taxonomic Levels"
# output_files_dict["third_page_header"] = "3.2 Taksonomik Seviyelere Göre Bulunma Yüzdeleri "


output_files_dict["third_page_text"] = ("With the abundance analysis, the bacteria contained in the sample were examined at different taxonomic levels. "
                                       "The abundances of bacteria at each taxonomic level in the figures below are shown on the pie chart.")
# output_files_dict["third_page_text"] = ("Bulunma yüzdelerine göre numunede bulunan bakteriler farklı taksonomik seviyelerde incelenmiştir. Aşağıdaki şekillerde her bir taksonomik seviyedeki bakteri bollukları pasta grafik üzerinde gösterilmiştir.")


create_pdf_inputs(output_files_dict, OUTPUT_PDF)
