# Grid refinement simple app

library(shiny)
library(ggplot2)
library(data.table)

# read in example datasets
simple_data <<- data.frame("X" = c(5, 4, 6, 2),
                           "Y" = c(10, 3, 7, 3))
#read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vT5dXhUhLHM93c23pbfnq9k2OApr6ZfV5Gjv1PklX_woOHkSZSDJ_kqMDfaoyxFm4Z5CzltOBeNvO_p/pub?gid=1459563270&single=true&output=csv")

complex_data <<- data.frame("X" = runif(n = 5, min = 1000, max = 10000),
                            "Y" = runif(n = 5, min = 1000, max = 10000))
  #read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQIJDP3KqUtMA3xsMrNxgXrNicpX1Pj6lJK3HSdyWbIRNwU2vImSf4DZTA3QV1ArqSFFpYgRox0NLdX/pub?gid=1285731657&single=true&output=csv")
Point_ID <<- 0

# Grid Functions
Make_Grid <- function(x_coords, y_coords, cell_size, buffer){
  #check for same length in x and y 
  if (length(x_coords) != length(y_coords)){
    return('The x and y coordinates do not have the same number of values')
  }
  # calculate the origin
  x_origin = floor(min(x_coords)) - (cell_size * buffer)
  y_origin = floor(min(y_coords)) - (cell_size * buffer)
  # calculate the farthest distance out
  x_dist = ceiling(max(x_coords)) + (cell_size * buffer)
  y_dist = ceiling(max(y_coords)) + (cell_size * buffer)
  # grid check, if last value in the sequence is not greater than distance in x and y need to add one more cell size to vector
  x_seq <- seq(x_origin ,x_dist, cell_size)
  if (max(x_seq) <= (x_dist - cell_size/2)){
    x_seq <- c(x_seq, max(x_seq) + cell_size)
  }
  y_seq <- seq(y_origin ,y_dist, cell_size)
  if (max(y_seq) <= (y_dist - cell_size/2)){
    y_seq <- c(y_seq, max(y_seq) + cell_size)
  }
  # create grid data set
  grid_data <- CJ(x = x_seq, y = y_seq)
  # add cell information
  grid_data$cell_size = cell_size
  grid_data$cell_id = 1:nrow(grid_data)
  grid_data
}

Ghost_Nodes <- function(Grid, x_coords, y_coords, ref_method){
  # method selection
  if (ref_method == 'Single') {
    Ghosts <- lapply(1:length(x_coords), function(a){
      temp_point_x <- x_coords[a]
      temp_point_y <- y_coords[a]
      temp_cell_size <- min(Grid$cell_size)
      # locate cell ids
      Grid_Nodes <- Grid[between(x = x, lower = temp_point_x-temp_cell_size/2, upper = temp_point_x+temp_cell_size/2, incbounds = T)][between(x = y, lower = temp_point_y-temp_cell_size/2, upper = temp_point_y+temp_cell_size/2, incbounds = T)]
      
      Grid_Nodes <- Grid_Nodes[cell_size == min(Grid_Nodes$cell_size)][1]
      
      Grid_Nodes
    })
  }else{
    
    if(ref_method != 'Radial'){return("No method found with that name. Please check spelling of the reference method")}
    # run for each point pair
    Ghosts <- lapply(1:length(x_coords), function(a){
      temp_point_x <- x_coords[a]
      temp_point_y <- y_coords[a]
      temp_cell_size <- min(Grid$cell_size)
      
      #find max/min x and y by multiplying the point location by 1.5 the current minimum distance in the grid
      min_x = temp_point_x - temp_cell_size*1.5
      min_y = temp_point_y - temp_cell_size*1.5
      max_x = temp_point_x + temp_cell_size*1.5
      max_y = temp_point_y + temp_cell_size*1.5
      
      # extract cell ids that will be needed
      Grid_Nodes <- Grid[between(x, min_x, max_x)][between(y, min_y, max_y)]
      Grid_Nodes
    })
  }
  Grid_Nodes <- rbindlist(Ghosts)
  Grid_Nodes
}

Quad_Tree <- function(Grid, ghost_nodes){
  # this is the outer loop that controls the number of refinements to do
  # using ghost nodes refine grid down number of steps desired
  new_cells <- lapply(ghost_nodes$cell_id, function(b) {
    temp_grid <- Grid[cell_id == b]
    # reduce each side by 2
    new_cell1 = data.table(x = temp_grid$x + temp_grid$cell_size/4, y = temp_grid$y + temp_grid$cell_size/4, cell_size = temp_grid$cell_size/2, cell_id = 0)
    new_cell2 = data.table(x = temp_grid$x + temp_grid$cell_size/4, y = temp_grid$y - temp_grid$cell_size/4, cell_size = temp_grid$cell_size/2, cell_id = 0)
    new_cell3 = data.table(x = temp_grid$x - temp_grid$cell_size/4, y = temp_grid$y - temp_grid$cell_size/4, cell_size = temp_grid$cell_size/2, cell_id = 0)
    new_cell4 = data.table(x = temp_grid$x - temp_grid$cell_size/4, y = temp_grid$y + temp_grid$cell_size/4, cell_size = temp_grid$cell_size/2, cell_id = 0)
    new_cells <- rbindlist(list(new_cell1,new_cell2, new_cell3, new_cell4))
  })
  # remove old cells and add in new cells, order by x and then y
  new_cells <- rbindlist(new_cells)
  Grid <- Grid[!cell_id %in% ghost_nodes$cell_id]
  Grid <- rbind(Grid, new_cells)
  Grid <- Grid[order(x,y)]
  Grid$cell_id <- 1:nrow(Grid)
  Grid
}

Grid_Refinement_Function <- function(x_coords, y_coords, cell_size, buffer, num_ref, ref_method) {
  temp_grid <- Make_Grid(x_coords = x_coords, y_coords = y_coords, cell_size = cell_size, buffer = buffer)
  # check if refinement number is zero and return normal grid
  if (num_ref > 0){
    # now we enter into the refine grid loop
    for (i in 1:num_ref){
      # first find the ghost nodes
      ghosts <- Ghost_Nodes(Grid = temp_grid, x_coords = x_coords, y_coords = y_coords, ref_method = ref_method)
      # then run quad tree refinement
      temp_grid <- Quad_Tree(Grid = temp_grid, ghost_nodes = ghosts)
    }
  }
  # if no more refinements return the updated grid
  temp_grid
}

Base_R_Lines <- function(grid_value) {
  
  # find all coordinates
  x1_3 <- grid_value$x - grid_value$cell_size/2
  x2_4 <- grid_value$x + grid_value$cell_size/2
  
  y1_2 <- grid_value$y + grid_value$cell_size/2
  y3_4 <- grid_value$y - grid_value$cell_size/2
  
  # make rect vectors
  x_vect <- c(x1_3,x2_4,x2_4,x1_3,x1_3)
  y_vect <- c(y1_2,y1_2,y3_4,y3_4,y1_2)
  
  # return as dataframe
  return(data.frame(x = x_vect, y = y_vect))
}


ui <- fluidPage(
  
  # App title ----
  titlePanel("Grid Refinement App"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      #fileInput(inputId = 'point_file', label = 'Select CSV File to Load', 
                #accept = c("text/csv", "text/comma-separated-values, .csv")),
      selectInput(inputId = "example_data", label = "Select Example Dataset", 
                  choices = c("Simple", "Complex")),
      #selectInput(inputId = 'x_coords', label = 'Choose X Coordinate Column', choices = NULL),
      #helpText('Be sure to select the correct column to use'),
      #selectInput(inputId = 'y_coords', label = 'Choose Y Coordinate Column', choices = NULL),
      #helpText('Be sure to select the correct column to use'),
      #numericInput(inputId = 'x_point', label = 'X Value', value = 0),
      #numericInput(inputId = 'y_point', label = 'Y Value', value = 0),
      numericInput(inputId = 'cell_size', label = 'Grid Size (Cell Size)', value = 100, min = 0),
      numericInput(inputId = 'buffer', label = 'Buffer', value = 2, min = 0),
      numericInput(inputId = 'num_ref', label = 'Number of Grid Refinments', value = 0, min = 0),
      selectInput(inputId = 'ref_method', label = 'Grid Refinement Method', choices = c('Single', 'Radial')),
      helpText('Single refines the grid at the piont location/Radial refines at the point and the surrounding cells'),
      #actionButton(inputId = 'add_button', label = 'Add Points'),
      actionButton(inputId = 'run_refinement', label = 'Create Grid'),
      #actionButton(inputId = 'help_button', label = 'Help')
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: grid graph
      plotOutput(outputId = 'Point_Graph',height = 1000, width = '100%')
      
    )
  )
)

# server section
server <- function(input, output, session) {
  
  # load example datasets
  observeEvent(eventExpr = input$example_data, 
               handlerExpr = {
                 if (input$example_data == "Simple") {
                   table <- simple_data
                   table <- cbind(data.table('Point ID' = 1:nrow(table)), table)
                   Point_Table <<- table
                 }
                 
                 if (input$example_data == "Complex") {
                   table <- complex_data
                   table <- cbind(data.table('Point ID' = 1:nrow(table)), table)
                   Point_Table <<- table
                 }
                 
               })
  
  # Help Button Pressed
  observeEvent(eventExpr = input$help_button, {
    showModal(modalDialog(
      includeHTML(path = 'Help_Text.html'),
      easyClose = T
    ))
  })
  
  # Add Button Pressed
  observeEvent(eventExpr = input$add_button, {
    # Error message if they try to add points to existing file
    if (nrow(Point_Table) > 1) {
      showModal(modalDialog(
        "You can not add points to an existing dataset at this time",
        easyClose = T
      ))
      return()
    }
    # increase ID number
    Point_ID = Point_ID + 1
    # check if table is created yet, if not create, if so bind new values
    if(Point_ID == 1){
      if (length(Point_Table) != 0){
        Point_Table <<- rbind(Point_Table, data.table('Point ID' = Point_ID, 'X ' = input$x_point, 'Y' = input$y_point))
      }else{
        Point_Table <<- data.table('Point ID' = Point_ID, 'X' = input$x_point, 'Y' = input$y_point)
      }
    }else{
      Point_Table <<- rbind(Point_Table, data.table('Point ID' = Point_ID, 'X' = input$x_point, 'Y' = input$y_point))
    }
    # change the x and y columns to the correct variable name
    updateSelectInput(session, inputId = 'x_coords', choices = colnames(Point_Table), selected = 'X')
    updateSelectInput(session, inputId = 'y_coords', choices = colnames(Point_Table), selected = 'Y')
    # plot table
    output$Point_Table <- DT::renderDataTable(expr = DT::datatable(data = Point_Table))
  })
  
  # Observe if table is active, if so then produce table
  observeEvent(eventExpr = input$point_file, ignoreNULL = T, ignoreInit = T, {
    #check if columns have label X and Y
    table <- input$point_file
    table <- fread(input = table$datapath)
    table <- cbind(data.table('Point ID' = 1:nrow(table)), table)
    output$Point_Table <- DT::renderDataTable(expr = DT::datatable(data = table))
    Point_Table <<- table
    
    #name_check <- colnames(table)
    #if (any(name_check == 'X') & any(name_check == 'Y')){
    # create table if file is uploaded
    #table <- cbind(data.table('Point ID' = 1:nrow(table)), table)
    #output$Point_Table <- DT::renderDataTable(expr = DT::datatable(data = table))
    #Point_Table <- table
    #}else{
    #validate("Please Upload Table with column names X and Y")
    #}
    
  })
  
  # read in table columns and update column names in inputs
  observe({
    updateSelectInput(session, inputId = 'x_coords', choices = colnames(Point_Table))
    updateSelectInput(session, inputId = 'y_coords', choices = colnames(Point_Table))
  })
  
  # If create grid button is clicked produce the graph and grid
  observeEvent(eventExpr = input$run_refinement, ignoreInit = T, ignoreNULL = T,{
    withProgress(message = "Creating Grid", value = 0, {
      incProgress(0.2)
      table_data = as.data.frame(Point_Table)
      x_coord = table_data[,"X"]
      y_coord = table_data[,"Y"]
      incProgress(0.2)
      # get grid
      grid_data <<- Grid_Refinement_Function(x_coords =  x_coord, y_coords = y_coord, 
                                             cell_size = input$cell_size, 
                                             buffer = input$buffer, 
                                             num_ref = input$num_ref, ref_method = input$ref_method)
    
      incProgress(0.2)
      # make base R rect dataset
      baseR_rects <- lapply(1:nrow(grid_data), FUN = function(x) Base_R_Lines(grid_value = grid_data[x]))
  
      # make new table just for points
      point_table <<- data.frame(x_coord,y_coord)
      incProgress(0.2)
    
      # render graph
      output$Point_Graph <- renderPlot(expr = {
        plot(point_table, 
             xlim = c(min(grid_data$x), max(grid_data$x)), 
             ylim = c(min(grid_data$y), max(grid_data$y)),
             col = "steelblue",
             pch = 16,
             cex = 2
             )
        lapply(baseR_rects, FUN = function(x) {
          
          lines(x$x, x$y)
        })
        incProgress(0.2)
      })
    })
  }) 
}

# Run the application 
shinyApp(ui = ui, server = server)
