test_that("weightedAdjMatrix returns erros when expected", {
    expect_error(weightedAdjMatrix("dog", "kegg", outputDir = "data/test.rds"),"Requested species or database currently not supported by `grahpite`. Run `pathwayDatabses`
        to get databases available.")
    expect_error(weightedAdjMatrix("hsapiens", "kegg", outputDir = "data/test.rds", beta = c(0,1,2)), "Beta has wrong length or names. See details for requirements")
    expect_error(weightedAdjMatrix("hsapiens", "kegg", outputDir = "data/test.rds", pathwayName = c("A", "B")), "Pathway names provided not detected in retrieved database")

    })

test_that("weightedAdjMatrix writes file to designated file path", {
    outputDir <- paste(tempdir(), "/test.rda", sep = "")
    expect_false(file.exists(outputDir))
    weightedAdjMatrix(species = "hsapiens",
                      database = "kegg",
                      outputDir = outputDir,
                      pathwayName = c("Glycolysis / Gluconeogenesis","Citrate cycle (TCA cycle)","Pentose phosphate pathway" ))
    expect_true(file.exists(outputDir))
    load(outputDir)
    expect_true(exists("BminsI"))
    expect_identical(names(BminsI), c("Glycolysis / Gluconeogenesis","Citrate cycle (TCA cycle)","Pentose phosphate pathway"))
    expect_true(stringr::str_detect(rownames(BminsI[[1]])[1], "ENTREZID:"))
    unlink(outputDir)

})
