library(signifinder)
library(testthat)
suppressPackageStartupMessages(library(SummarizedExperiment))
set.seed(12345)

test_that("interferonSign based on Barkley's work", {
    rmatrix <- .fakeData("Interferon_Barkley")
    malign <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
    myres <- interferonSign(rmatrix, isMalignant = malign)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Interferon_Barkley" %in% colnames(colData(myres)))
    expect_length(colData(myres)[, "Interferon_Barkley"], ncol(assay(myres)))
    expect_type(colData(myres)[, "Interferon_Barkley"], "double")
    expect_message(interferonSign(rmatrix, isMalignant = malign), "100")
})

test_that("EMTSign based on Miow's work", {
    pyrnames <- c("EMT_Miow_Epithelial", "EMT_Miow_Mesenchymal")
    pname <- sample(pyrnames, 1)
    rmatrix <- .fakeData(pname)
    expect_warning(myres <- EMTSign(
        rmatrix, normalize = FALSE), "less than 30% of its genes")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(pname %in% colnames(colData(myres)))
    expect_length(colData(myres)[, pname], ncol(assay(myres)))
    expect_type(colData(myres)[, pname], "double")
    expect_message(expect_warning(EMTSign(rmatrix, normalize = FALSE)), "100")
})

test_that("EMTSign based on Mak's work", {
    rmatrix <- .fakeData("EMT_Mak")
    myresMak <- EMTSign(rmatrix, author = "Mak")
    expect_true(is(myresMak, "SummarizedExperiment"))
    expect_true("EMT_Mak" %in% colnames(colData(myresMak)))
    expect_length(colData(myresMak)[, "EMT_Mak"], ncol(assay(myresMak)))
    expect_type(colData(myresMak)[, "EMT_Mak"], "double")
    expect_message(EMTSign(rmatrix, author = "Mak"), "100")
})

test_that("EMTSign based on Cheng's work", {
    rmatrix <- .fakeData("EMT_Cheng")
    myres <- EMTSign(rmatrix, author = "Cheng")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("EMT_Cheng" %in% colnames(colData(myres)))
    expect_length(colData(myres)[, "EMT_Cheng"], ncol(assay(myres)))
    expect_type(colData(myres)[, "EMT_Cheng"], "double")
    expect_message(EMTSign(rmatrix, author = "Cheng"), "100")
})

test_that("EMTSign based on Thompson's work", {
  rmatrix <- .fakeData("EMT_Thompson")
  myres <- EMTSign(rmatrix, author = "Thompson")
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("EMT_Thompson" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "EMT_Thompson"], ncol(assay(myres)))
  expect_type(colData(myres)[, "EMT_Thompson"], "double")
  expect_message(EMTSign(rmatrix, author = "Thompson"), "100")
})

test_that("EMTSign based on Barkley's work", {
  rmatrix <- data.frame()
  PanBarkley <- SignatureNames[grep("Barkley", SignatureNames)]
  for(i in PanBarkley){
    rmatrix <- rbind(rmatrix, as.data.frame(.fakeData(i)))}
  malign <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
  myres <- EMTSign(rmatrix, inputType = "sc", author = "Barkley",
                   isMalignant = malign)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("EMT_Barkley_cEMT" %in% colnames(colData(myres)))
  expect_length(colData(myres)[,"EMT_Barkley_pEMT"],ncol(assay(myres)))
  expect_type(colData(myres)[,"EMT_Barkley_cEMT"], "double")
  expect_message(EMTSign(rmatrix, inputType = "sc", author = "Barkley",
                         isMalignant = malign), "100")
})

test_that("hypoxiaSign based on Barkley's work", {
  rmatrix <- data.frame()
  PanBarkley <- SignatureNames[grep("Barkley", SignatureNames)]
  for(i in PanBarkley){
    rmatrix <- rbind(rmatrix, as.data.frame(.fakeData(i)))}
  malign <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
  myres <- hypoxiaSign(rmatrix, inputType = "sc", author = "Barkley",
                   isMalignant = malign)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("Hypoxia_Barkley" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "Hypoxia_Barkley"], ncol(assay(myres)))
  expect_type(colData(myres)[, "Hypoxia_Barkley"], "double")
  expect_message(hypoxiaSign(rmatrix, inputType = "sc", author = "Barkley",
                             isMalignant = malign), "100")
})

test_that("ASCSign work", {
    rmatrix <- .fakeData("ASC_Smith")
    myres <- ASCSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("ASC_Smith" %in% colnames(colData(myres)))
    expect_length(colData(myres)[, "ASC_Smith"], ncol(assay(myres)))
    expect_type(colData(myres)[, "ASC_Smith"], "double")
    expect_message(ASCSign(rmatrix), "100")
})

test_that("chemokineSign work", {
    rmatrix <- .fakeData("Chemokines_Messina")
    myres <- chemokineSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("Chemokines_Messina" %in% colnames(colData(myres)))
    expect_length(colData(myres)[, "Chemokines_Messina"], ncol(assay(myres)))
    expect_type(colData(myres)[, "Chemokines_Messina"], "double")
    expect_message(chemokineSign(rmatrix), "100")
})

test_that("PassONSign work", {
    rmatrix <- .fakeData("PassON_Du")
    myres <- PassONSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("PassON_Du" %in% colnames(colData(myres)))
    expect_length(colData(myres)[, "PassON_Du"], ncol(assay(myres)))
    expect_type(colData(myres)[, "PassON_Du"], "double")
})

test_that("CISSign work", {
    rmatrix <- .fakeData("CIS_Robertson")
    myres <- CISSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("CIS_Robertson" %in% colnames(colData(myres)))
    expect_length(colData(myres)[, "CIS_Robertson"], ncol(assay(myres)))
    expect_type(colData(myres)[, "CIS_Robertson"], "double")
    expect_message(CISSign(rmatrix), "100")
})

test_that("HRDSSign work", {
    rmatrix <- .fakeData("HRDS_Lu")
    myres <- HRDSSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("HRDS_Lu" %in% colnames(colData(myres)))
    expect_length(colData(myres)[, "HRDS_Lu"], ncol(assay(myres)))
    expect_type(colData(myres)[, "HRDS_Lu"], "double")
    expect_message(HRDSSign(rmatrix), "100")
})

test_that("DNArepSign work", {
    rmatrix <- .fakeData("DNArep_Kang")
    myres <- DNArepSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("DNArep_Kang" %in% colnames(colData(myres)))
    expect_length(colData(myres)[, "DNArep_Kang"], ncol(assay(myres)))
    expect_type(colData(myres)[, "DNArep_Kang"], "double")
    expect_message(DNArepSign(rmatrix), "100")
})

test_that("IPRESSign work", {
    rmatrix <- .fakeData("IPRES_Hugo")
    myres <- IPRESSign(rmatrix)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("IPRES_Hugo" %in% colnames(colData(myres)))
    expect_length(colData(myres)[, "IPRES_Hugo"], ncol(assay(myres)))
    expect_type(colData(myres)[, "IPRES_Hugo"], "double")
})

test_that("ECM work", {
    pyrnames <- c("ECM_Chakravarthy_up", "ECM_Chakravarthy_down")
    pname <- sample(pyrnames, 1)
    rmatrix <- .fakeData(pname)
    expect_warning(myres <- ECMSign(
        rmatrix, normalize = FALSE), "less than 30% of its genes")
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true(pname %in% colnames(colData(myres)))
    expect_length(colData(myres)[, pname], ncol(assay(myres)))
    expect_type(colData(myres)[, pname], "double")
    expect_message(expect_warning(ECMSign(rmatrix, normalize = FALSE)), "100")
})

test_that("IPSOVSign work", {
    rmatrix  <- .fakeData("IPSOV_Shen")
    myres <- IPSOVSign(rmatrix, normalize = FALSE)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("IPSOV_Shen" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"IPSOV_Shen"], ncol(assay(myres)))
    expect_type(colData(myres)[,"IPSOV_Shen"], "double")
    expect_message(IPSOVSign(rmatrix, normalize = FALSE), "100")
})

test_that("stateSign based on Neftel's work", {
    rmatrix <- data.frame()
    for(i in SignatureNames){
        rmatrix <- rbind(rmatrix, as.data.frame(.fakeData(i)))}
    malign <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
    myres <- stateSign(rmatrix, author = "Neftel",
                       isMalignant = malign)
    expect_true(is(myres, "SummarizedExperiment"))
    expect_true("State_Neftel_OPC" %in% colnames(colData(myres)))
    expect_length(colData(myres)[,"State_Neftel_AC"],ncol(assay(myres)))
    expect_type(colData(myres)[,"State_Neftel_MES1"], "double")
    expect_message(stateSign(rmatrix, author = "Neftel",
                             isMalignant = malign), "100")
})

test_that("stateSign based on Barkley's work", {
  rmatrix <- data.frame()
  PanBarkley <- SignatureNames[grep("Barkley", SignatureNames)]
  for(i in PanBarkley){
    rmatrix <- rbind(rmatrix, as.data.frame(.fakeData(i)))}
  malign <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
  myres <- stateSign(rmatrix, author = "Barkley",
                   isMalignant = malign)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("State_Barkley_AC" %in% colnames(colData(myres)))
  expect_length(colData(myres)[,"State_Barkley_OPC"],ncol(assay(myres)))
  expect_type(colData(myres)[,"State_Barkley_NPC"], "double")
  expect_message(stateSign(rmatrix, author = "Barkley",
                           isMalignant = malign), "100")
})

test_that("TinflamSign based on Thompson's work", {
  rmatrix <- .fakeData("Tinflam_Thompson")
  myres <- TinflamSign(rmatrix, author = "Thompson")
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("Tinflam_Thompson" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "Tinflam_Thompson"], ncol(assay(myres)))
  expect_type(colData(myres)[, "Tinflam_Thompson"], "double")
  expect_message(TinflamSign(rmatrix, author = "Thompson"), "100")
})

test_that("cellCycleSign based on Barkley's work", {
  rmatrix <- .fakeData("CellCycle_Barkley")
  malign <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
  myres <- cellCycleSign(rmatrix, inputType = "sc", author = "Barkley",
                       isMalignant = malign)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("CellCycle_Barkley" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "CellCycle_Barkley"], ncol(assay(myres)))
  expect_type(colData(myres)[, "CellCycle_Barkley"], "double")
  expect_message(cellCycleSign(rmatrix, inputType = "sc", author = "Barkley",
                             isMalignant = malign), "100")
})

test_that("CombinedSign work", {
  rmatrix <- .fakeData("Combined_Thompson")
  myres <- CombinedSign(rmatrix)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("Combined_Thompson" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "Combined_Thompson"], ncol(assay(myres)))
  expect_type(colData(myres)[, "Combined_Thompson"], "double")
  expect_message(CombinedSign(rmatrix), "100")
})

test_that("stateSign based on Tirosh's work", {
  rmatrix <- data.frame()
  for(i in SignatureNames){
    rmatrix <- rbind(rmatrix, as.data.frame(.fakeData(i)))}
  malign <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
  myres <- stateSign(rmatrix, author = "Tirosh",
                     isMalignant = malign)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("State_Tirosh_MITF" %in% colnames(colData(myres)))
  expect_length(colData(myres)[,"State_Tirosh_AXL"],ncol(assay(myres)))
  expect_type(colData(myres)[,"State_Tirosh_AXL"], "double")
  expect_message(stateSign(rmatrix, author = "Tirosh",
                           isMalignant = malign), "100")
})

test_that("APMSign based on Thompson's work", {
  rmatrix <- .fakeData("APM_Thompson")
  myres <- APMSign(rmatrix, author = "Thompson")
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("APM_Thompson" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "APM_Thompson"], ncol(assay(myres)))
  expect_type(colData(myres)[, "APM_Thompson"], "double")
  expect_message(APMSign(rmatrix, author = "Thompson"), "100")
})

test_that("APMSign based on Wang's work", {
  rmatrix <- .fakeData("APM_Wang")
  myres <- APMSign(rmatrix)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("APM_Wang" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "APM_Wang"], ncol(assay(myres)))
  expect_type(colData(myres)[, "APM_Wang"], "double")
  expect_message(APMSign(rmatrix), "100")
})

test_that("MPSSign work", {
  rmatrix <- .fakeData("MPS_PerezGuijarro")
  myres <- MPSSign(rmatrix)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("MPS_PerezGuijarro" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "MPS_PerezGuijarro"], ncol(assay(myres)))
  expect_type(colData(myres)[, "MPS_PerezGuijarro"], "double")
  expect_message(MPSSign(rmatrix), "100")
})

test_that("IRGSign work", {
  rmatrix <- .fakeData("IRG_Yang")
  myres <- IRGSign(rmatrix)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("IRG_Yang" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "IRG_Yang"], ncol(assay(myres)))
  expect_type(colData(myres)[, "IRG_Yang"], "double")
  expect_message(IRGSign(rmatrix), "100")
})

test_that("TGFBSign work", {
  rmatrix <- .fakeData("TGFB_Mariathasan")
  myres <- TGFBSign(rmatrix)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("TGFB_Mariathasan" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "TGFB_Mariathasan"], ncol(assay(myres)))
  expect_type(colData(myres)[, "TGFB_Mariathasan"], "double")
  expect_message(TGFBSign(rmatrix), "100")
})

test_that("ADOSign work", {
  rmatrix <- .fakeData("ADO_Sidders")
  myres <- ADOSign(rmatrix)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("ADO_Sidders" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "ADO_Sidders"], ncol(assay(myres)))
  expect_type(colData(myres)[, "ADO_Sidders"], "double")
  expect_message(ADOSign(rmatrix), "100")
})

test_that("MITFlowPTENnegSign work", {
  rmatrix <- .fakeData("MITFlowPTENneg_Cabrita")
  myres <- MITFlowPTENnegSign(rmatrix)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("MITFlowPTENneg_Cabrita" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "MITFlowPTENneg_Cabrita"], ncol(assay(myres)))
  expect_type(colData(myres)[, "MITFlowPTENneg_Cabrita"], "double")
  expect_message(MITFlowPTENnegSign(rmatrix), "100")
})

test_that("LRRC15CAFSign work", {
  rmatrix <- .fakeData("LRRC15CAF_Dominguez")
  myres <- LRRC15CAFSign(rmatrix)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("LRRC15CAF_Dominguez" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "LRRC15CAF_Dominguez"], ncol(assay(myres)))
  expect_type(colData(myres)[, "LRRC15CAF_Dominguez"], "double")
  expect_message(LRRC15CAFSign(rmatrix), "100")
})

test_that("SCSubtypeSign work", {
  pyrnames <- c(
      "SCSubtype_Wu_Basal", "SCSubtype_Wu_Her2E",
      "SCSubtype_Wu_LumA", "SCSubtype_Wu_LumB")
  pname <- sample(pyrnames, 1)
  rmatrix <- .fakeData(pname)
  malign <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
  expect_warning(expect_warning(expect_warning(myres <- SCSubtypeSign(
      rmatrix, isMalignant = malign), "less than 30% of its genes")))
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true(pname %in% colnames(colData(myres)))
  expect_length(colData(myres)[,pname],ncol(assay(myres)))
  expect_type(colData(myres)[,pname], "double")
  expect_message(expect_warning(expect_warning(expect_warning(
      SCSubtypeSign(rmatrix, isMalignant = malign)))), "100")
})

test_that("ICBResponseSign work", {
  pyrnames <- c("ICBResponse_Chen_responder", "ICBResponse_Chen_nonresponder")
  pname <- sample(pyrnames, 1)
  rmatrix <- .fakeData(pname)
  expect_warning(expect_warning(
      myres <- ICBResponseSign(rmatrix), "less than 30% of its genes"))
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true(pname %in% colnames(colData(myres)))
  expect_length(colData(myres)[, pname], ncol(assay(myres)))
  expect_type(colData(myres)[, pname], "double")
  expect_message(expect_warning(expect_warning(ICBResponseSign(rmatrix))), "100")
})

test_that("COXISSign work", {
  rmatrix <- .fakeData("COXIS_Bonavita")
  myres <- COXISSign(rmatrix)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("COXIS_Bonavita" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "COXIS_Bonavita"], ncol(assay(myres)))
  expect_type(colData(myres)[, "COXIS_Bonavita"], "double")
  expect_message(COXISSign(rmatrix), "100")
})

test_that("stressSign based on Barkley's work", {
  rmatrix <- .fakeData("Stress_Barkley")
  malign <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
  myres <- stressSign(rmatrix, isMalignant = malign)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("Stress_Barkley" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "Stress_Barkley"], ncol(assay(myres)))
  expect_type(colData(myres)[, "Stress_Barkley"], "double")
  expect_message(stressSign(rmatrix, isMalignant = malign), "100")
})

test_that("oxphosSign based on Barkley's work", {
  rmatrix <- data.frame()
  PanBarkley <- SignatureNames[grep("Barkley", SignatureNames)]
  for(i in PanBarkley){
    rmatrix <- rbind(rmatrix, as.data.frame(.fakeData(i)))}
  malign <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
  myres <- oxphosSign(rmatrix, isMalignant = malign)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("Oxphos_Barkley" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "Oxphos_Barkley"], ncol(assay(myres)))
  expect_type(colData(myres)[, "Oxphos_Barkley"], "double")
  expect_message(oxphosSign(rmatrix, isMalignant = malign), "100")
})

test_that("metalSign based on Barkley's work", {
  rmatrix <- data.frame()
  PanBarkley <- SignatureNames[grep("Barkley", SignatureNames)]
  for(i in PanBarkley){
    rmatrix <- rbind(rmatrix, as.data.frame(.fakeData(i)))}
  malign <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
  myres <- metalSign(rmatrix, isMalignant = malign)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("Metal_Barkley" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "Metal_Barkley"], ncol(assay(myres)))
  expect_type(colData(myres)[, "Metal_Barkley"], "double")
  expect_message(metalSign(rmatrix, isMalignant = malign), "100")
})

test_that("CD39CD8TcellSign work", {
  rmatrix <- .fakeData("CD39CD8Tcell_Chow")
  myres <- CD39CD8TcellSign(rmatrix)
  expect_true(is(myres, "SummarizedExperiment"))
  expect_true("CD39CD8Tcell_Chow" %in% colnames(colData(myres)))
  expect_length(colData(myres)[, "CD39CD8Tcell_Chow"], ncol(assay(myres)))
  expect_type(colData(myres)[, "CD39CD8Tcell_Chow"], "double")
  expect_message(CD39CD8TcellSign(rmatrix), "100")
})
