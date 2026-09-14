#' Evaluate Downstream Model Performance After Imputation
#' 
#' Generic function to evaluate downstream model performance after DIMV imputation.
#' Model-agnostic: supports any model that accepts a formula and data frame, 
#' or a pre-trained model object.
#' 
#' @param imputed_data Data frame results from dimv_impute_new() or dimv_impute_multiple()
#' @param formula Formula object, e.g. y ~ x1 + x2
#' @param model_method Character string specifying the modeling method:
#'   - "xgboost": Use xgboost package (requires xgboost installed)
#'   - "ranger": Use ranger's random forest implementation (optional)
#'   - "randomForest": Use randomForest (optional)
#'   - "glmnet": Use elastic-net regression (optional)
#'   - "glm": Use generalized linear model (base R)
#'   - "lm": Use linear model (base R, default)
#'   - "custom": Use a pre-specified model object
#' @param model Optional pre-trained model object (used when model_method = "custom")
#' @param test_idx Optional vector of row indices for test data; if NULL, uses 80/20 split
#' @param ... Additional arguments passed to the modeling function
#' 
#' @return A list containing:
#'   \item{model}{The trained model object}
#'   \item{predictions}{Predicted values on test data}
#'   \item{actual}{Actual values on test data}
#'   \item{mse}{Mean squared error on test data}
#'   \item{method}{Character string of the method used}
#'   \item{formula}{The formula used for modeling}
#' 
#' @export
#' @importFrom stats as.formula glm lm model.matrix
#' @examples
#' library(dimvR)
#' set.seed(123)
#' X <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
#' for (j in seq_along(X)) X[sample(1:nrow(X), 10), j] <- NA
#' imp <- dimv_train(X, lambda = 0.1)
#' X_imp <- dimv_impute_new(imp, X)
#' X_imp$y <- rnorm(nrow(X_imp))
#' 
#' # Example with lm (base R, no extra packages needed)
#' result <- evaluate_downstream(X_imp, y ~ x1 + x2, model_method = "lm")
#' result\$mse
#' 
#' # Example with xgboost (if available)
#' \dontrun{
#' result <- evaluate_downstream(X_imp, y ~ x1 + x2, model_method = "xgboost")
#' result\$mse
#' }
evaluate_downstream <- function(imputed_data, formula, 
                                model_method = c("lm", "glm", "ranger", "randomForest",
                                                  "glmnet", "xgboost", "custom"),
                                model = NULL, test_idx = NULL, ...) {
  
  model_method <- match.arg(model_method)
  dots <- list(...)

  if (!is.data.frame(imputed_data)) {
    imputed_data <- as.data.frame(imputed_data)
  }
  
  # Split data into train/test if test_idx not provided
  if (is.null(test_idx)) {
    n <- nrow(imputed_data)
    set.seed(123)  # For reproducible split
    test_idx <- sample(seq_len(n), size = ceiling(0.2 * n))
  }
  
  test_data <- imputed_data[test_idx, , drop = FALSE]
  train_data <- imputed_data[-test_idx, , drop = FALSE]
  
  # Handle formula supplied as text and preserve its terms for train/test.
  f <- formula
  if (inherits(f, "character")) {
    f <- as.formula(f)
  }
  
  # Get the response variable name
  response_name <- all.vars(f)[1]
  
  train_formula <- f
  train_x <- stats::model.matrix(train_formula, data = train_data)[, -1, drop = FALSE]
  test_x <- stats::model.matrix(train_formula, data = test_data)[, -1, drop = FALSE]
  train_y <- train_data[[response_name]]
  
  # Get method
  method <- model_method
  
  # Train model based on method
  model <- switch(method,
    
    # Linear model (base R, always available)
    "lm" = {
      lm(train_formula, data = train_data)
    },
    
    # GLM (base R)
    "glm" = {
      glm(train_formula, data = train_data, ...)
    },

    # ranger (optional, fast random forest implementation)
    "ranger" = {
      if (!requireNamespace("ranger", quietly = TRUE)) {
        stop("Package 'ranger' needed for this method. Install it or choose another method.")
      }
      ranger::ranger(formula = train_formula, data = train_data,
                     respect.unordered.factors = "order", ...)
    },

    # randomForest (optional, established random forest implementation)
    "randomForest" = {
      if (!requireNamespace("randomForest", quietly = TRUE)) {
        stop("Package 'randomForest' needed for this method. Install it or choose another method.")
      }
      randomForest::randomForest(formula = train_formula, data = train_data, ...)
    },

    # Elastic-net regression (optional)
    "glmnet" = {
      if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package 'glmnet' needed for this method. Install it or choose another method.")
      }
      # `s` is a prediction-time choice in glmnet. Remove it from fitting
      # arguments if supplied, then use it below for one prediction column.
      fit_args <- dots
      fit_args$s <- NULL
      do.call(glmnet::glmnet, c(list(x = train_x, y = train_y), fit_args))
    },
    
    # XGBoost (requires package)
    "xgboost" = {
      # Check if xgboost is available
      if (!requireNamespace("xgboost", quietly = TRUE)) {
        stop("Package 'xgboost' needed for this method. Install it or choose another method.")
      }
      
      # Prepare xgboost matrix
      dtrain <- xgboost::xgboost(data = as.matrix(train_x),
                                   label = train_y,
                                   nrounds = 100,
                                   objective = "reg:squarederror",
                                   ...)
    },
    
    # Custom model (pre-trained)
    "custom" = {
      if (is.null(model)) {
        stop("Must provide a pre-trained model object when method = 'custom'")
      }
      model
    }
  )
  
  test_y <- test_data[[response_name]]
  
  # Make predictions
  predictions <- switch(method,
    "lm" = predict(model, newdata = test_data),
    "glm" = predict(model, newdata = test_data),
    "ranger" = {
      as.numeric(stats::predict(model, data = test_data)$predictions)
    },
    "randomForest" = {
      as.numeric(stats::predict(model, newdata = test_data))
    },
    "glmnet" = {
      s_value <- dots$s
      if (is.null(s_value)) s_value <- model$lambda[1L]
      as.numeric(stats::predict(model, newx = test_x, s = s_value))
    },
    "xgboost" = {
      # For xgboost, we need the model matrix
      if (!requireNamespace("xgboost", quietly = TRUE)) {
        stop("Package 'xgboost' needed for xgboost method")
      }
      predict(model, newdata = as.matrix(test_x))
    },
    "custom" = predict(model, newdata = test_data)
  )
  
  # Calculate MSE
  mse <- mean((predictions - test_y)^2, na.rm = TRUE)
  
  list(
    model = model,
    predictions = predictions,
    actual = test_y,
    mse = mse,
    method = method,
    formula = as.character(formula),
    test_idx = test_idx
  )
}
