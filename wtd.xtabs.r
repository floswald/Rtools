

# todo: add weights argument in first line
# introduce if (is.null(weights)) switch
# if null, current function
# else, find a way to extract and carry around weights

wtd.xtabs <- function (formula = ~., data = parent.frame(), weights, sparse = FALSE, 
    na.action, exclude = c(NA, NaN), drop.unused.levels = FALSE) 
{
    if (missing(formula) && missing(data)) 
        stop("must supply either 'formula' or 'data'")
    if (!missing(formula)) {
        formula <- as.formula(formula)
        if (!inherits(formula, "formula")) 
            stop("'formula' missing or incorrect")
    }
    if (any(attr(terms(formula, data = data), "order") > 1)) 
        stop("interactions are not allowed")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m$... <- m$exclude <- m$drop.unused.levels <- m$sparse <- NULL
    m[[1L]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    if (length(formula) == 2L) {
        by <- mf
        y <- NULL
    }
    else {
        i <- attr(attr(mf, "terms"), "response")
        by <- mf[-i]
        y <- mf[[i]]
    }
    by <- lapply(by, function(u) {
        if (!is.factor(u)) 
            u <- factor(u, exclude = exclude)
        u[, drop = drop.unused.levels]
    })
    if (!sparse) {
        x <- if (is.null(y)) 
            do.call("table", by)
        else if (NCOL(y) == 1L) 
            tapply(y, by, sum)
        else {
            z <- lapply(as.data.frame(y), tapply, by, sum)
            array(unlist(z), dim = c(dim(z[[1L]]), length(z)), 
                dimnames = c(dimnames(z[[1L]]), list(names(z))))
        }
        x[is.na(x)] <- 0
        class(x) <- c("xtabs", "table")
        attr(x, "call") <- match.call()
        x
    }
    else {
        if (length(by) != 2L) 
            stop("xtabs(*, sparse=TRUE) applies only to two-way tables")
        if (is.null(tryCatch(loadNamespace("Matrix"), error = function(e) NULL))) 
            stop("xtabs(*, sparse=TRUE) needs package \"Matrix\" correctly installed")
        if (length(i.ex <- unique(unlist(lapply(by, function(f) which(is.na(f))))))) 
            by <- lapply(by, `[`, -i.ex)
        rows <- by[[1L]]
        cols <- by[[2L]]
        rl <- levels(rows)
        cl <- levels(cols)
        if (is.null(y)) 
            y <- rep.int(1, length(rows))
        as(new("dgTMatrix", i = as.integer(rows) - 1L, j = as.integer(cols) - 
            1L, x = as.double(y), Dim = c(length(rl), length(cl)), 
            Dimnames = list(rl, cl)), "CsparseMatrix")
    }
}

