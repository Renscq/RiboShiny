# install.packages("shinydashboard")
# library(shiny)
# library(shinydashboard)



library(shiny)
library(shinydashboard)
library(httpuv)

# UI部分
ui <- dashboardPage(
  dashboardHeader(title = "Client IP Address Dashboard"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard"))
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "dashboard",
              fluidRow(
                box(title = "Client Information", status = "primary", solidHeader = TRUE, 
                    width = 12,
                    textOutput("ip_address")  # 显示IP地址
                )
              )
      )
    )
  )
)

server <- function(input, output, session) {
  # browser()
  
  client_ip <- reactive({
    browser()
    
    ip <- tryCatch({
      session$request$HTTP_HOST
      session$clientData$url_hostname
      
      # paste0(session$clientData$url_protocol, 
      #        session$clientData$url_pathname, session$clientData$url_pathname, 
      #        session$clientData$url_hostname, ':', session$clientData$url_port)
      
    }, error = function(e) {
      "Unknown IP"
    })
    ip
  })
  
  # show the ip address
  output$ip_address <- renderText({
    paste("Your IP address is:", client_ip())
  })
}

# 运行Shiny应用
shinyApp(ui, server)
