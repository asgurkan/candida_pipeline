import html_sections as hs

def create_html_report(html_dict, output_file):
    
    html_string = """
    <!DOCTYPE html>
    <html>

    <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css" type="text/css">
    <link rel="stylesheet" href="styles/theme.css" type="text/css">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    </head>

    <body style="box-shadow: 0px 0px 40px  black;">
    <nav class="navbar navbar-expand-lg bg-dark navbar-dark sticky-top" style="">
        <div class="container"> <button class="navbar-toggler navbar-toggler-right border-0 p-0" type="button" data-toggle="collapse" data-target="#navbar14">
            <p class="navbar-brand mb-0 text-white">
            <i class="fa d-inline fa-lg fa-stop-circle"></i> Massive Bioinformatics </p>
        </button>
        <div class="collapse navbar-collapse justify-content-start" id="navbar14">
            <ul class="navbar-nav ml-auto">
            <li class="nav-item mx-1"> <a class="nav-link" target="_blank" href="https://www.massivebioinformatics.com">
                <img src="styles/logoMassive.png" height="50" width="50">
                </a> </li>
            </ul>
        </div>
        </div>
    </nav>
    <div class="py-5 text-center" style="	background-image: url(styles/bacteria.jpg);	background-size: cover;	background-position: top left;	background-repeat: repeat;">
        <div class="container-fluid">
        <div class="row">
            <div class="bg-light text-dark p-5 mx-auto col-md-8 col-10 shadow" style="	box-shadow: 0px 0px 40px  black; opacity:0.6">
            <h2 class="display-3 text-dark" style="	text-shadow: 0px 0px 0px black;">Comparative Metagenomic Analysis Report</h2>
            <p class="mb-3 lead">Massive Bioinformatics<br>info@massivebioinformatics.com</p>
            </div>
        </div>
        </div>
    </div>
    <div class="py-3 text-center" id="main_content">
        <div class="container">
        <p style="page-break-before: always" ></p>
        """ + hs.test_infos(html_dict['metadata'], part=1)+ """
        <p style="page-break-before: always" ></p>
        
        <div class="py-3 text-center">
        <div class="container">
        """ + hs.bioinformatic(part=2)+ """
        """ + hs.bioinformatic_krona(html_dict)+ """
        <p style="page-break-before: always" ></p>
        """ + hs.bioinformatic_genus(html_dict)+ """
        <p style="page-break-before: always" ></p>
        """ + hs.bioinformatic_similarity_matrix(html_dict)+ """
        <p style="page-break-before: always" ></p>
        """ + hs.bioinformatic_alpha(html_dict)+ """
        <p style="page-break-before: always" ></p>
        """ + hs.bioinformatic_pca(html_dict)+ """
        <p style="page-break-before: always" ></p>
        """ + hs.bioinformatic_beta(html_dict)+ """
        <p style="page-break-before: always" ></p>
        """ + hs.bioinformatic_lefse()+ """
        <p style="page-break-before: always" ></p>
        </div>
        </div>
        """ + hs.test_metot(part=3)+ """
        <p style="page-break-before: always" ></p>
        """ + hs.raw_data(html_dict, part=4)+ """
        <p style="page-break-before: always" ></p>


        
    </div>    
    </div>
    </body>
    </html>    
    """

    print(output_file.split(".")[0]+".html")
    f = open(output_file.split(".")[0]+".html",'w')
    f.write(html_string)
    f.close()