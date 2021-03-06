context("Test: getDatasets()")


equal_lists <- function(df1, df2) {
        isTRUE(all.equal(df1, df2))
}


test_that("The getDatasets() interface works properly..",{
        
        skip_on_cran()
        skip_on_travis()
        
    expect_output(marts <- getMarts()[ , "mart"])
    expect_output(getDatasets(mart = "ENSEMBL_MART_ENSEMBL"))
    
})





