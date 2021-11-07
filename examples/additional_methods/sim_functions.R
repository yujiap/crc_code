# library(umap)     # UMAP
# library(e1071)    # KNN
# library(ranger)   # RF
# library(snifter)  # tSNE. Requires R 4.1.0
# library(kernlab)  # SVM
# library(class)    # kNN
# library(caret)

# NOTES:
# UMAP uses nearest neighbor embedding
# tNSE as additional columns: https://mark-borg.github.io/blog/2016/tsne-ml/
# Interpolation based t-SNE for >2 dimensions is currently unsupported (and generally a bad idea)

train_test_LDA = function(Ztrn, Ttrn, Ztst, Ttst, fctr_tst=NULL) {
  if (ncol(Ztrn) > nrow(Ztrn)) stop("Number of features exceeds number of observations.")
  fit = lda(Ztrn, Ttrn)
  pred = predict(fit, Ztst)$class
  acc = mean(pred == Ttst)
  subgroup_acc = if (!is.null(fctr_tst)) subgroup_acc(pred, Ttst, fctr_tst)
  return(list(acc = acc, subgroup_acc = subgroup_acc))
}

train_test_SVM = function(Ztrn, Ttrn, Ztst, Ttst, params, fctr_tst=NULL) {
  method = params$method
  colnames(Ztrn) = colnames(Ztst) = paste0("V", 1:ncol(Ztrn))
  train_control = trainControl(method = "cv", number = 10, returnData = FALSE, trim = TRUE, returnResamp = "none")
  fit = if (method == "radial") {
    train(x = Ztrn, y = Ttrn, method = "svmRadial", trControl = train_control, preProcess = c("center","scale"), tuneLength = 10)
  } else if (method == "poly") {
    train(x = Ztrn, y = Ttrn, method = "svmPoly", trControl = train_control, preProcess = c("center","scale"), tuneLength = 5) # for computational speed
  } else if (method == "linear") {
    train(x = Ztrn, y = Ttrn, method = "svmLinear", trControl = train_control, preProcess = c("center","scale"), tuneGrid = expand.grid(C = seq(0, 2, length = 9)))
  }
  pred = predict.train(fit, Ztst)
  acc = mean(pred == Ttst)
  subgroup_acc = if (!is.null(fctr_tst)) subgroup_acc(pred, Ttst, fctr_tst)
  return(list(acc = acc, subgroup_acc = subgroup_acc))
}

train_test_tsneLDA = function(Ztrn, Ttrn, Ztst, Ttst, params=NULL, fctr_tst=NULL) {
  Ztrnproj = fitsne(Ztrn, n_components = 2L)
  Ztstproj = project(x = Ztrnproj, new = Ztst, old = Ztrn)
  # Objects should be of class matrix (first), not snifter
  Ztrnproj = matrix(Ztrnproj, ncol = 2)
  Ztstproj = matrix(Ztstproj, ncol = 2)
  return(train_test_LDA(Ztrnproj, Ttrn, Ztstproj, Ttst, fctr_tst))
}

train_test_tsneSVM = function(Ztrn, Ttrn, Ztst, Ttst, params, fctr_tst=NULL) {
  Ztrnproj = fitsne(Ztrn, n_components = 2L)
  Ztstproj = project(x = Ztrnproj, new = Ztst, old = Ztrn)
  # Objects should be of class matrix (first), not snifter
  Ztrnproj = matrix(Ztrnproj, ncol = 2)
  Ztstproj = matrix(Ztstproj, ncol = 2)
  out = train_test_SVM(Ztrnproj, Ttrn, Ztstproj, Ttst, params = params, fctr_tst = fctr_tst)
  return(list(acc = out$acc, subgroup_acc = out$subgroup_acc))
}

train_test_tsnekNN = function(Ztrn, Ttrn, Ztst, Ttst, params, fctr_tst=NULL) {
  Ztrnproj = fitsne(Ztrn, n_components = 2L)
  Ztstproj = project(x = Ztrnproj, new = Ztst, old = Ztrn)
  # Objects should be of class matrix (first), not snifter
  Ztrnproj = matrix(Ztrnproj, ncol = 2)
  Ztstproj = matrix(Ztstproj, ncol = 2)
  out = train_test_fastkNN(Ztrnproj, Ttrn, Ztstproj, Ttst, fctr_tst = fctr_tst)
  return(list(acc = out$acc, subgroup_acc = out$subgroup_acc))
}

train_test_umapLDA = function(Ztrn, Ttrn, Ztst, Ttst, params=NULL, fctr_tst=NULL) {
  umapObj = umap(Ztrn)
  Ztrnproj = umapObj$layout
  Ztstproj = predict(umapObj, Ztst)
  return(train_test_LDA(Ztrnproj, Ttrn, Ztstproj, Ttst, fctr_tst))
}

train_test_fastkNN = function(Ztrn, Ttrn, Ztst, Ttst, params=NULL, fctr_tst=NULL) {
  svdZtrn = svd(Ztrn)
  Ztrn1 = svdZtrn$u %*% diag(svdZtrn$d)
  Ztst1 = Ztst %*% svdZtrn$v
  colnames(Ztrn1) = colnames(Ztst1) = paste0("V", 1:ncol(Ztrn1))
  train_control = trainControl(method = "cv", number = 10, returnData = FALSE, trim = TRUE, returnResamp = "none")
  fit = train(x = Ztrn1, y = Ttrn, method = "knn", trControl = train_control, tuneLength = 10) # preProcess = c("center", "scale")
  pred = predict.train(fit, Ztst1)
  acc = mean(pred == Ttst)
  subgroup_acc = if (!is.null(fctr_tst)) subgroup_acc(pred, Ttst, fctr_tst)
  return(list(acc = acc, subgroup_acc = subgroup_acc))
}

train_test_umapSVM = function(Ztrn, Ttrn, Ztst, Ttst, params, fctr_tst=NULL) {
  umapObj = umap(Ztrn)
  Ztrnproj = umapObj$layout
  Ztstproj = predict(umapObj, Ztst)
  out = train_test_SVM(Ztrnproj, Ttrn, Ztstproj, Ttst, params = params, fctr_tst = fctr_tst)
  return(list(acc = out$acc, subgroup_acc = out$subgroup_acc))
}

train_test_umapkNN = function(Ztrn, Ttrn, Ztst, Ttst, params, fctr_tst=NULL) {
  umapObj = umap(Ztrn)
  Ztrnproj = umapObj$layout
  Ztstproj = predict(umapObj, Ztst)
  out = train_test_fastkNN(Ztrnproj, Ttrn, Ztstproj, Ttst, fctr_tst = fctr_tst)
  return(list(acc = out$acc, subgroup_acc = out$subgroup_acc))
}

train_test_RF = function(Ztrn, Ttrn, Ztst, Ttst, params, fctr_tst=NULL) {
  fit = ranger(x = Ztrn, y = Ttrn, classification = TRUE)
  pred = predict(fit, Ztst)
  acc = mean(pred$predictions==Ttst)
  subgroup_acc = if (!is.null(fctr_tst)) subgroup_acc(pred, Ttst, fctr_tst)
  return(list(acc = acc, subgroup_acc = subgroup_acc))
}

train_test_fsvaDlda = function(Ztrn, Ttrn, Ztst, Ttst, params=NULL, fctr_tst=NULL) {
  if (sum(apply(Ztrn, 2, sd) == 0) > 0) stop("One or more feature variances are constant.")
  trn_mod = model.matrix(~factor(Ttrn))
  trn_mod0 = rep(1, length(Ttrn))
  trn_sv = sva(t(Ztrn), trn_mod, trn_mod0)
  fsvaobj = fsva(t(Ztrn), trn_mod, trn_sv, t(Ztst))
  fit = myDlda(t(fsvaobj$db), Ttrn, prior = c(0.5, 0.5))
  out = train_test_myDlda(Ztrn, Ttrn, Ztst, Ttst, params=NULL, fctr_tst=NULL)
  return(list(acc = out$acc, subgroup_acc = out$subgroup_acc))
}
