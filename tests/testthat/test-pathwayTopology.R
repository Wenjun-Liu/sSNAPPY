test_that("retrieve_topology returns erros when expected", {
    expect_error(retrieve_topology(database = "random"),"Requested database currently not supported. See example for how to find databases available")
    expect_error(retrieve_topology(database = "kegg", beta = c(0,1,2)), "Beta has wrong length or names. See details for requirements")
    expect_error(retrieve_topology(database = "kegg", pathwayName = c("A", "B")), "Pathway names provided not detected in retrieved database")

    })

test_that("retrieve_topology returns expected output", {
    temp <- retrieve_topology(database = "kegg",
                      pathwayName = c("Glycolysis / Gluconeogenesis","Citrate cycle (TCA cycle)","Pentose phosphate pathway" ))
    expect_identical(names(temp), c("Glycolysis / Gluconeogenesis","Citrate cycle (TCA cycle)","Pentose phosphate pathway"))
    expect_true(stringr::str_detect(rownames(temp[[1]])[1], "ENTREZID:"))

})
