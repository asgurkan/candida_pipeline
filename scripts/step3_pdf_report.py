from fpdf import FPDF
from PIL import Image

WIDTH = 210
HEIGHT = 297

def get_image(img):
    im = Image.open(img)
    return im


TEXT_SIZE = 10
HEADER1_TEXT_SIZE = 15
HEADER2_TEXT_SIZE = 13
HEADER1_PARAPH = 13
HEADER2_PARAPH = 15
TEXT_PARAPH = 15
TABLE_COLUMN_TEXT_SIZE = 11
TABLE_VALUE_TEXT_SIZE = 11
PAGE_BACKGROUND = "gorseller/page_eng.png"

class PDF(FPDF):
    def header(self):
        self.image("gorseller/page.png", 0, 0, WIDTH, HEIGHT)

    def footer(self):
        self.set_y(-15)
        self.set_font("helvetica", "I", 12)
        self.cell(0, 15, f"{self.page_no()}", align="R")
    
def call_text_font(pdf, given_size, given_style):
    pdf.add_font("PTSans", "", "gorseller/PTSans/PTSans-Regular.ttf", uni=True)
    pdf.add_font("PTSans", "B", "gorseller/PTSans/PTSans-Bold.ttf", uni=True)
    pdf.set_font("PTSans", given_style, size=given_size)
    pdf.set_text_color(64, 64, 64) 

def render_table_header(pdf, col_width, line_height, TABLE_COL_NAMES):
    call_text_font(pdf, TABLE_COLUMN_TEXT_SIZE, "B")
    for col_name in TABLE_COL_NAMES:
        pdf.cell(col_width, line_height, col_name, align="C", border=1)
    pdf.ln(line_height)
    pdf.set_font(style="")  # disabling bold text
    call_text_font(pdf, TABLE_VALUE_TEXT_SIZE, "")


def render_table_header_for_sixth_page(pdf, col_widths, line_height, TABLE_COL_NAMES):
    call_text_font(pdf, TABLE_COLUMN_TEXT_SIZE, "B")
    for i, col_name in enumerate(TABLE_COL_NAMES):
        pdf.cell(col_widths[i], line_height, col_name, align="C", border=1)
    pdf.ln(line_height)
    pdf.set_font(style="")  # disabling bold text
    call_text_font(pdf, TABLE_VALUE_TEXT_SIZE, "")


def first_page(first_page_dict, pdf):
    pdf.add_page()
    pdf.image("gorseller/page1.png", 0, 0, WIDTH, HEIGHT)
    # pdf.image("gorseller/avitek.png", 150, 250)
    call_text_font(pdf, 16,"")
    pdf.set_xy(42,217)


def second_page(second_page_dict,pdf):
    pdf.add_page()
    pdf.image(PAGE_BACKGROUND,0,0,WIDTH, HEIGHT)
    call_text_font(pdf, HEADER1_TEXT_SIZE, "B")
    pdf.set_xy(13,50)
    pdf.write(15, second_page_dict["header"])
    call_text_font(pdf, TEXT_SIZE,"")
    pdf.set_xy(TEXT_PARAPH, 65)
    pdf.multi_cell(w=180, h=6, txt=second_page_dict["text1"], markdown=True)
    pdf.set_xy(15,110)
    pdf.multi_cell(w=180, h=6, txt=second_page_dict["text2"], markdown=True)


def quality_page(quality_page_dict, pdf):
    pdf.add_page()
    pdf.image(PAGE_BACKGROUND, 0, 0, WIDTH, HEIGHT)

    call_text_font(pdf, 13, "B")
    pdf.set_xy(HEADER1_PARAPH, 50)
    pdf.write(13, quality_page_dict["header"])

    call_text_font(pdf, 10, "")
    pdf.set_xy(TEXT_PARAPH, 65)
    line_height = pdf.font_size * 2.5

    pdf.multi_cell(w=180, h=6, txt=quality_page_dict["text1"], markdown=True)
    pdf.ln()

    quality_statistics = quality_page_dict["quality_statistics"]
    line_height = pdf.font_size * 2.5
    col_width = pdf.epw / len(quality_statistics.columns)
    pdf.set_draw_color(172, 208, 255)

    render_table_header(pdf, col_width, line_height, quality_statistics.columns)

    for row_i, row_v in quality_statistics.iterrows():
        if pdf.will_page_break(line_height):
            pdf.add_page()
            render_table_header(pdf, col_width, line_height, quality_statistics.columns)
        for item in row_v:
            pdf.cell(col_width, line_height, str(item), border=1, align="C")
        pdf.ln(line_height)

    pdf.image(get_image(quality_page_dict["quality_per_seq_image_path"]), 50, 150, 120, 70)
    pdf.set_xy(80,230)
    call_text_font(pdf, 10, "B")
    pdf.write(4, "Figure 1. Quality per-Seq Graph")


def alpha_diversity_page(alpha_diversity_page_dict, pdf):
    pdf.add_page()
    pdf.image(PAGE_BACKGROUND, 0, 0, WIDTH, HEIGHT)

    call_text_font(pdf, HEADER1_TEXT_SIZE, "B")
    pdf.set_xy(HEADER1_PARAPH, 50)
    pdf.write(15, alpha_diversity_page_dict["header1"])

    call_text_font(pdf, HEADER2_TEXT_SIZE, "B")
    pdf.set_xy(HEADER2_PARAPH, 65)
    pdf.write(13, alpha_diversity_page_dict["header2"])

    call_text_font(pdf, TEXT_SIZE, "")
    pdf.set_xy(TEXT_PARAPH, 75)
    line_height = pdf.font_size * 2.5

    pdf.multi_cell(w=180, h=6, txt=alpha_diversity_page_dict["text1"], markdown=True)
    pdf.ln()

    alpha_diversity_table = alpha_diversity_page_dict["alpha_diversity_table"]
    line_height = pdf.font_size * 2.5
    col_width = pdf.epw / len(alpha_diversity_table.columns)
    pdf.set_draw_color(172, 208, 255)
    render_table_header(pdf, col_width, line_height, alpha_diversity_table.columns)
    
    for row_i, row_v in alpha_diversity_table.iterrows():
        if pdf.will_page_break(line_height):
            pdf.add_page()
            render_table_header(pdf, col_width, line_height, alpha_diversity_table.columns)
        for item in row_v:
            pdf.cell(col_width, line_height, str(item), border=1, align="C")
    
    pdf.ln()
    call_text_font(pdf, TEXT_SIZE, "")
    pdf.set_xy(TEXT_PARAPH, 120)
    pdf.multi_cell(w=180, h=6, txt=alpha_diversity_page_dict["text2"], markdown=True)
    

def third_page(third_page_dict, pdf): #5
    pdf.add_page()
    pdf.image(PAGE_BACKGROUND, 0, 0, WIDTH, HEIGHT)

    call_text_font(pdf, HEADER2_TEXT_SIZE, "B")
    pdf.set_xy(HEADER2_PARAPH, 50)
    pdf.write(13, third_page_dict["header"])

    call_text_font(pdf, TEXT_SIZE,"")
    pdf.set_xy(TEXT_PARAPH,65)
    pdf.multi_cell(w=180, h=6, txt=third_page_dict["text"], markdown=True)

    call_text_font(pdf, TEXT_SIZE,"")
    pdf.image(get_image(third_page_dict["phylum"]), 60, 90, 120, 70)
    pdf.set_xy(80,165)
    call_text_font(pdf, 10, "B")
    pdf.write(4, "Figure 2. Pie Chart in Phylum level")
    # pdf.write(4, "Figür 2. Şube seviyesinde pasta grafiği")
    pdf.image(get_image(third_page_dict["class"]), 60, 175, 120, 70)
    pdf.set_xy(80,250)
    call_text_font(pdf, 10, "B")
    pdf.write(4, "Figure 3. Pie Chart in Class level")
    # pdf.write(4, "Figür 3. Sınıf seviyesinde pasta grafiği")


def forth_page(forth_page_dict, pdf): #6
    pdf.add_page()
    pdf.image(PAGE_BACKGROUND, 0, 0, WIDTH, HEIGHT)

    call_text_font(pdf, TEXT_SIZE,"")
    pdf.image(get_image(forth_page_dict["order"]), 55, 60, 120, 70)
    pdf.set_xy(80,140)
    call_text_font(pdf, 10, "B")
    pdf.write(4, "Figure 4. Pie Chart in Order level")
    # pdf.write(4, "Figür 4. Takım seviyesinde pasta grafiği")
    pdf.image(get_image(forth_page_dict["family"]), 55, 160, 120, 70)
    pdf.set_xy(80,240)
    call_text_font(pdf, 10, "B")
    pdf.write(4, "Figure 5. Pie Chart in Family level")
    # pdf.write(4, "Figür 5. Familya seviyesinde pasta grafiği")

def fifth_page(fifth_page_dict, pdf): #7
    pdf.add_page()
    pdf.image(PAGE_BACKGROUND, 0, 0, WIDTH, HEIGHT)

    call_text_font(pdf, TEXT_SIZE,"")
    pdf.image(get_image(fifth_page_dict["genus"]), 46, 60, 120, 70)
    pdf.set_xy(80, 140)
    call_text_font(pdf, 10, "B")
    pdf.write(4, "Figure 6. Pie Chart in Genus level")
    # pdf.write(4, "Figür 6. Cins seviyesinde pasta grafiği")


    pdf.image(get_image(fifth_page_dict["species"]), 55, 160, 120, 70)
    pdf.set_xy(80, 240)
    call_text_font(pdf, 10, "B")
    pdf.write(4, "Figure 7. Pie Chart in Species level")
    # pdf.write(4, "Figür 7. Tür seviyesinde pasta grafiği")


def sixth_page(sixth_page_dict, pdf): #8
    pdf.add_page()
    pdf.image(PAGE_BACKGROUND, 0, 0, WIDTH, HEIGHT)

    call_text_font(pdf, 13, "B")
    pdf.set_xy(13, 50)
    pdf.write(13, "3.3 ABUNDANCE TABLE" )
    # pdf.write(13, "3.3 BULUNMA YÜZDELERİ TABLOSU" )

    call_text_font(pdf, 10, "")
    pdf.set_xy(10, 65)
    line_height = pdf.font_size * 2.5
    
    # Define column widths here
    col_widths = [100, 28, 28, 28] 
    pdf.set_draw_color(172, 208, 255)

    render_table_header_for_sixth_page(pdf, col_widths, line_height, sixth_page_dict.columns)
    
    for row_i, row_v in sixth_page_dict.iterrows():
        if pdf.will_page_break(line_height):
            pdf.add_page()
            pdf.set_xy(10, 50)
            render_table_header_for_sixth_page(pdf, col_widths, line_height, sixth_page_dict.columns)

        taxon = str(row_v["Taxon"])
        taxonomy = str(row_v["Taxonomy"])
        read_count = str(row_v["Read Count"])
        abundance = str(row_v["Abundance (%)"])

        # taxon = str(row_v["Takson"])
        # taxonomy = str(row_v["Taksonomi"])
        # read_count = str(row_v["Okuma Sayısı"])
        # abundance = str(row_v["Bolluk (%)"])

        pdf.cell(col_widths[0], line_height, taxon, align="L", border=1)
        pdf.cell(col_widths[1], line_height, taxonomy, align="C", border=1)
        pdf.cell(col_widths[2], line_height, read_count, align="C", border=1)
        pdf.cell(col_widths[3], line_height, abundance, align="C", border=1) 
        pdf.ln(line_height)


def create_report(report_content, output_pdf):
    pdf = PDF()
    first_page(report_content["sample_name"], pdf)
    second_page(report_content["second_page"], pdf)
    quality_page(report_content["quality_page"], pdf)
    alpha_diversity_page(report_content["alpha_diversity_page"], pdf)
    third_page(report_content["third_page"], pdf)
    forth_page(report_content["forth_page"], pdf)
    fifth_page(report_content["fifth_page"], pdf)
    sixth_page(report_content["raw_data"], pdf)
    pdf.output(output_pdf)