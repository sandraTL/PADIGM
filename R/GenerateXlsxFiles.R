
#assoDataXlsx <- function(x){

# exporting data.frame to excel is easy with xlsx package
sheetname <- "mysheet"
write.xlsx(realData1, "mydata.xlsx", sheetName=sheetname)
file <- "mydata.xlsx"
# but we want to highlight cells if value greater than or equal to 5
wb <- loadWorkbook(file)              # load workbook
fo1 <- Fill(foregroundColor="blue")   # create fill object # 1
cs1 <- CellStyle(wb, fill=fo1)        # create cell style # 1
fo2 <- Fill(foregroundColor="red")    # create fill object # 2
cs2 <- CellStyle(wb, fill=fo2)        # create cell style # 2
sheets <- getSheets(wb)               # get all sheets
sheet <- sheets[[sheetname]]          # get specific sheet
rows <- getRows(sheet, rowIndex=2:(nrow(realData1)+1))     # get rows

# 1st row is headers

cells <- getCells(rows, colIndex=1:3)         # get cells
print(cells)
# in the wb I import with loadWorkbook, numeric data starts in column 2
# The first column is row numbers.  The last column is "yes" and "no" entries, so
# we do not include them, thus we use colIndex = 2:cols

values <- lapply(cells, getCellValue) # extract the cell values

highlightblue <- NULL
for (i in names(values)) {

    x <- as.numeric(values[i])
    if (x>5 & !is.na(x)) {
        print(i)
        highlightblue <- c(highlightblue, i)
    }
}

# find cells meeting conditional criteria < 5
highlightred <- NULL
for (i in names(values)) {
    x <- as.numeric(values[i])
    if (x<5 & !is.na(x)) {
        print(i)
        highlightred <- c(highlightred, i)
    }
}

 lapply(names(cells[highlightblue]),
        function(ii)setCellStyle(cells[[ii]],cs1))

 lapply(names(cells[highlightred]),
        function(ii)setCellStyle(cells[[ii]],cs2))

saveWorkbook(wb, file)
#}
