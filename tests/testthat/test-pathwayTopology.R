test_that("retrieve_topology returns erros when expected", {
    expect_error(retrieve_topology(database = "random"))
    expect_error(retrieve_topology(database = "kegg", species = "random"))
    expect_error(retrieve_topology(database = "kegg", species = "hsapien", beta = c(0,1,2)), "Beta has wrong length or names. See details for requirements")
    expect_error(retrieve_topology(database = "kegg", species = "hsapien", keyword = c("random")))
})

test_that("retrieve_topology returns expected output", {
    temp <- retrieve_topology(database = "kegg", species = "hsapien",
                              keyword = "estrogen")
    expect_equal(names(temp), "kegg.Estrogen signaling pathway")
    expect_true(stringr::str_detect(rownames(temp[[1]])[1], "ENTREZID:"))

    # if multiple dataset were provided
    temp <- retrieve_topology(database = c("kegg", "path"), species = "hsapien",
                              keyword = "estrogen")
    expect_equal(names(temp), c("kegg.Estrogen signaling pathway","pathbank.Androgen and Estrogen Metabolism"))
    expect_true(stringr::str_detect(rownames(temp[[1]])[1], "ENTREZID:"))
})
