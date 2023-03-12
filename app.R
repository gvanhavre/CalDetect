#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(sf)
library(sp)
library(ggplot2)
library(ggforce)
library(leaflet)
library(lwgeom)

# Carregar Poligonos
cnsa <- st_read("http://portal.iphan.gov.br/geoserver/ows?service=wfs&version=2.0.0&request=GetFeature&typename=SICG%3Asitios_pol&srsName=EPSG%3A4674")
cnsa <- st_cast(cnsa, "GEOMETRYCOLLECTION") %>%  st_collection_extract("POLYGON")
sf_use_s2(FALSE)

# Calcular diametro
cnsa$area <- st_area(cnsa)
cnsa$diametro <- 2*sqrt(cnsa$area/pi)

# Define UI for app´lication that draws plots
ui <- fluidPage(
  
  # Application title
  titlePanel(title=tags$div(img(src="title.png", width="100"), "Arq_CalDetect"),windowTitle ="Arq_CalDetect"),
  
  # Sidebar with a slider input
  sidebarLayout(
    sidebarPanel(
      h5("Ferramenta de auxílio à calibração do espaçamento entre sondagens em Arqueologia."),
      sliderInput("radius", "Raio de busca (km):",min=0,max=1000,value=100, step=10),
      sliderInput("user_lim", "Espaçamento (m):",min = 0,max = 500,value = 50, step=5),
      numericInput("user_lon", "Longitude (graus decimais):", 0, min=-180, max=180),
      numericInput("user_lat", "Latitude (graus decimais):", 0, min=-90, max=90),
      actionButton("do", "Calcular")
    ),
    
    mainPanel("",
              fluidRow(
                verticalLayout(
                  uiOutput("resumo"),
                  verbatimTextOutput("summary"),
                  uiOutput("resumo2"),
                  leafletOutput("mapa"), 
                  uiOutput("res"),
                  tableOutput("table"),
                  plotOutput("prob"),
                  uiOutput("ref")
              )
    )
  )
)
)

# Define server logic required to draw plots
server <- function(input, output) {
  observeEvent(input$do, {
    # Ponto de Origem
    temp <- data.frame(cbind(input$user_lon,input$user_lat))
    ptx <- st_as_sf(temp, coords=c("X1","X2"))
    st_crs(ptx) <- st_crs(cnsa)
    cnsa$dist <- st_distance(cnsa, ptx)
    buffer <- st_buffer(ptx, input$radius/111.120)
    
    # Preparar dados
    data_tmp <- subset(cnsa, as.numeric(cnsa$dist) < (input$radius*1000))
    
    # Média local
    if (nrow(data_tmp)>0) {
      for (i in 1:nrow(data_tmp)) {
        data_tmp$P[i] <- if (data_tmp$P[i] <- (pi*((as.numeric(data_tmp$diametro[i])^2))/(input$user_lim^2))*100 > 100) {100} else {(pi*((as.numeric(data_tmp$diametro[i])^2))/(input$user_lim^2))*100}
        data_tmp$count[i] <- 1
      }
    colnames(data_tmp)[18:20] <- c("Área", "Diâmetro", "Distância")
    x <- as.data.frame(table(cut(data_tmp$P,breaks=seq.int(from=min(data_tmp$P),to=max(data_tmp$P),by=1))))
    x$D <- seq(1,nrow(x),1)
    x2 <- as.data.frame(table(cut(data_tmp$P,breaks=10)))
    x2$Var1 <- seq(from=10, to=100, by=10)
    for (i in 1:nrow(x2)) {
      x2$Prop[i] <- round(x2$Freq[i]/sum(x2$Freq)*100)
    }
    colnames(x2) <- c("Intervalo", "Frequência", "Proporção")
    
    p3 <- ggplot(x, aes(x=D, y=Freq)) + 
      geom_bar(stat="identity", aes(fill = D)) + 
      scale_fill_gradientn(colours=c("darkred","white", "yelow", "white","darkgreen"), values=scales::rescale(c(min(data_tmp$P, na.rm=T),quantile(data_tmp$P)[2],median(data_tmp$P, na.rm=T),quantile(data_tmp$P)[4],max(data_tmp$P, na.rm=T)))) +
      xlab("Probabilidades por sítio na área de busca") +
      ylab("Frequência")
    }
    
    output$resumo <- renderUI(tags$div(
      tags$h3("Análise"),
      tags$h5("Coordenadas"), p(input$user_lon,"° Oeste e ",input$user_lat,"º Sul"),
      tags$h5("Espaçamento"), p(input$user_lim,"metros entre as sondagens."),
      tags$h5("Contexto"), p(input$radius,"quilômetros na volta do local, com um total de ", sum(data_tmp$count),"sítios."),
      tags$h4("Resumo dos dados")
    ))
    
    output$summary <- renderPrint({
      if (nrow(data_tmp)>0) {summary(data_tmp[,18:20])} else {NULL}
    })
    
    output$resumo2 <- renderUI(tags$div(
      tags$h4("Mapa de Localização")
    ))

    output$mapa <- renderLeaflet({
      leaflet(options = leafletOptions(minZoom = 4)) %>%
        addProviderTiles("OpenStreetMap") %>%
        setView( lng = input$user_lon
                 , lat = input$user_lat
                 , zoom = 4 ) %>%
        addCircles(lng=input$user_lon, lat=input$user_lat, radius = input$radius*1000, weight = 1, color = "darkred", fillColor = "red", fillOpacity = 0.1, popup="Area para calibração") %>%
        addCircles(lng=input$user_lon, lat=input$user_lat, radius = 100, weight = 1, color = "darkred", fillColor = "red", fillOpacity = 1)
    })
    
    output$res <- renderUI(tags$div(
      tags$h4("Resultados"),
      if (nrow(data_tmp)>0) {tags$div(
        p("Considerando as características dos",strong(sum(data_tmp$count)),strong("sítios"),"localizados numa",strong("área de"),strong(input$radius),strong("quilômetros"),"em torno das coordenadas inseridas, a metodologia proposta, com um",strong("espaçamento de"), strong(input$user_lim), strong("metros"),"entre sondagens, permite considerar que há uma",strong("probabilidade teórica de"), strong(round(median(data_tmp$P, na.rm=T),2)),strong("%"),"de encontrar, pelo menos, metade dos sítios semelhantes que apresentem as maiores dispersões de material."),
        p("A probabilidade representa apenas uma base para a avaliação da metodologia proposta. Está baseada na hipótese que os sítios conhecidos são semehantes àqueles que podem vir a ser encontrados. Novos contextos, com características diferentes, não entram em consideração. O resultado pode também variar de acordo com o raio de busca e a inclusão de novos sítios nas bases de dados."),
        p(strong("Em nenhum caso, este resultado pode servir de justificativa para não realizar estudos arqueológicos."))
        )}
      else {p("Sem dados. Tente novamente com um raio de busca maior.")}
    ))
    
    output$table <- renderTable(
      if (nrow(data_tmp)>0) {x2} else {NULL})
    
    output$prob <- renderPlot({
      if (nrow(data_tmp)>0) {p3} else {NULL}
    })
    
    output$ref <- renderUI(tags$div(
      tags$h4("Referências"),
      p("Nance, J. D. & Ball,B.F. 1986. No Surprises? The Reliability and Validity of Test Pit Sampling. In: American Antiquity, 51(3), p. 457-483."),
      p("Lightfoot, K. G. 1986. Regional Surveys in the Eastern United States: The Strengths and Weaknesses of Implementing Subsurface Testing Programs. In: American Antiquity, 51 (3), p. 484-504."),
      p("Wobst, H. M. 1983. We Can't See the Forest for the Trees: Sampling and the Shapes of Archaeological Distributions. In: Moore, J. A. & Keene, A. S. (Eds.) 1983. Archaeological Hammers and Theories. Londres: Academic Press, p. 37-85."),
      p("Araújo, A. 2001. Teoria e Método em Arqueologia Regional: Um Estudo de Caso no Alto Paranapanema, Estado de São Paulo. Tese de Doutorado. Museu de Arqueologia e Etnologia, Universidade de São Paulo."),
      p("IPHAN. 2023. SICG - Poligonais dos sítios arqueológicos. Brasília. Disponível em",tags$a("http://portal.iphan.gov.br/geoserver/ows?service=wfs&version=2.0.0&request=GetFeature&typename=SICG%3Asitios_pol&srsName=EPSG%3A4674"))
    ))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)