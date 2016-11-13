#require(sp)
#require(maptools)

acfTOshp <- function(file.GM      = NULL, 
                   name.SHP     = NULL,
                   exportPoints = T){
  GMapper <- read.table(file.GM, header = F, sep = ' ')
  grids <- grep(pattern = "DESCRIPTION=SCORE_", GMapper$V1)
  for (i in 1:length(grids)){
    grid <- as.data.frame(matrix(NA, 5, 2))
    grid[1, 1] <- as.numeric(as.character(GMapper[grids[i] + 7, ]))
    grid[1, 2] <- as.numeric(as.character(GMapper[grids[i] + 8, ]))
    grid[2, 1] <- as.numeric(as.character(GMapper[grids[i] + 9, ]))
    grid[2, 2] <- as.numeric(as.character(GMapper[grids[i] + 10, ]))
    grid[3, 1] <- as.numeric(as.character(GMapper[grids[i] + 11, ]))
    grid[3, 2] <- as.numeric(as.character(GMapper[grids[i] + 12, ]))
    grid[4, 1] <- as.numeric(as.character(GMapper[grids[i] + 13, ]))
    grid[4, 2] <- as.numeric(as.character(GMapper[grids[i] + 14, ]))
    grid[5, 1] <- as.numeric(as.character(GMapper[grids[i] + 15, ]))
    grid[5, 2] <- as.numeric(as.character(GMapper[grids[i] + 16, ]))
    colnames(grid) <- c('x', 'y')
    grid$y <- grid$y * (-1)
    txt.to.poly <- Polygon(grid)
    ID <- paste("ID", i, sep = "-")
    ID.poly <- Polygons(list(txt.to.poly), ID = ID)
    poly.to.class.SP <- SpatialPolygons(list(ID.poly))
    df <- data.frame(Area = i, row.names=ID)
    SP.df <- SpatialPolygonsDataFrame(poly.to.class.SP, df)
    proj4string(SP.df) <- "+proj=longlat +datum=WGS84"
    SP.df$Score <- sub('*DESCRIPTION=SCORE_*', '', GMapper[grids[i], ])
    RGB <- sub(' *BORDER_COLOR=RGB*', '', GMapper[grids[i] + 1, ])
    RGB <- gsub('[[:punct:]]', ' ', RGB)
    RGB <- unlist(strsplit(RGB, " "))
    SP.df$Red <- as.numeric(RGB[2])
    SP.df$Green <- as.numeric(RGB[3])
    SP.df$Blue<-as.numeric(RGB[4])
    writePolyShape(SP.df, paste(i, "polytemp", sep = "") )
    }
  shps <- list.files(pattern = "*polytemp.shp")
  poly0 <- NA
  #para cada shp una en un solo shp
  s <- '1polytemp.shp'
  for (s in shps) {
    #cargue el  poligono
    tem1 <- readShapePoly(s)
    poly1 <-  Polygons(tem1@polygons[[1]]@Polygons, ID = tem1@data$SP_ID)
    if (is.na(poly0)){
      #cargue el primer poligono
      poly.list <- poly1
      data1 <- tem1
      poly0 <- 1 
      }else{
        # este parte se hace para todos los poligono porque poly0 no es 1
        #sino cargue el segundo poligono
        # la info del segundo poligono asignela a poly2
        poly2 <- poly1
        #asigne el primer poligono a poly1
        poly1 <- poly.list
        # una el primer poligono con el segundo y asigneloa una lista
        poly.list <- as.list(c(poly1, poly2))
        # si poly es igual a 1, entonces es porque ya tiene la primera lista de 
        #poligonos (1 y 2)
        if (poly0 == 1) {
          # este parte se hace para solo para e?? segundo loop porque poly0 no es 1
          #asigne el segundo poligono a temp2
          tem2 <- tem1
          #asigne el primer poligono a temp1
          tem1 <- data1
          #ahora data1 sera la informacion de la lista de poligonos
          data1 <- rbind(tem1@data, tem2@data)
          #asine poly0 2 para que este paso no se repita
          poly0 <- 2
          }else{
            # este parte se hace para todos los poligono porque poly0 no es 1
            # si poly0 es diferente de 1 o 2
            # temp2 sera el siguiente poligono diferente alos poligonos
            #primero y al segundo
            tem2 <- tem1
            #temp1 sera lalista de poligonos, en la tercera ronda corresponde a la
            #lista del primero y el segundo poligono
            tem1 <- data1
            # este sera la lista de los poligonos juntos, en el tercer loop sera 
            # poligono 1, 2 y 3
            data1 <- rbind(tem1, tem2@data)
          }
      }
    }
  # cuando termine de unir en una lista todos los poligonos,pase a clase
  #SpatialPolygons 
  total.polys <- SpatialPolygons(poly.list)
  #organizar el orden de plot
  #total.polys@plotOrder<-total.polys@plotOrder[order(nchar(total.polys@plotOrder), total.polys@plotOrder)]
  # una la informacion de cada poligono
  rownames(data1) <- data1$SP_ID
  # pase a clase SpatialPolygonsDataFrame
  SP.df <- SpatialPolygonsDataFrame(total.polys, data1)
  #delete <- "rm *temp*"
  # borre todos los archivos
  unlink("*polytemp*", recursive = FALSE, force = FALSE)
  # system(delete)
  writePolyShape(SP.df, fn = name.SHP)
  # Escriba la table de las especies endemiscas de este poligono
  if (exportPoints==T){
    grids2 <- grep(pattern = "DESCRIPTION=", GMapper$V1)
    points <- grids2[which(!grids2 %in% grids)]
    table <- as.data.frame(matrix(data = NA,nrow = length(points), 3))
    for (j in 1:length(points)){
      colnames(table) <- c('species', 'decimalLongitude', 'decimalLatitude')
      table$species[j] <- sub('*DESCRIPTION=', '', GMapper[points[j], ])
      table$decimalLongitude[j] <- as.numeric(as.character(GMapper[points[j] + 3, ]))
      table$decimalLatitude[j] <- as.numeric(as.character(GMapper[points[j] + 4, ])) * (-1)
    }
    readAndWrite(action = 'write', frmt = 'saveTXT', path = NULL, 
               name = paste('EndemicSpecies', name.SHP, '.txt', sep='-'), object = table)
  
  }
}
