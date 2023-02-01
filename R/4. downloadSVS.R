
query.harmonized <- GDCquery(project = "TCGA-COAD", 
                             data.category = "Biospecimen", 
                             data.type = "Slide Image",
                             experimental.strategy = "Diagnostic Slide"
                             )  

query.harmonized  %>% 
  getResults %>% 
  head  %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))


query.harmonized <- GDCquery(project = "TCGA-COAD", 
                             data.category = "Biospecimen", 
                             data.type = "Slide Image",
                             experimental.strategy = "Diagnostic Slide",
                             data.format="SVS")  
GDCdownload(
  query = query.harmonized, 
  method = "api", 
  files.per.chunk = 10
)
