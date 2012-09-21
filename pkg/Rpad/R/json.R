json <- function (x, ...) 
  UseMethod("json")

json.list <- function(x, ...) {
  res <- array("", length(x))
  for (i in seq(len=length(x)))
    res[i] <- json(x[[i]])
  nms <- names(x)
  if (is.null(nms)) nms = seq(from = 0, len = length(x))
  idx <- which(nms == "")
  nms[idx] <- seq(from = 0, along = idx)
  result <- paste('{', paste('"', nms, '":', res, collapse = ',', sep = ''), '}', sep = '')
  class(result) <- "json"
  result
}

json.numeric <- function(x, ...) {
  if (!is.null(names(x)))
    result <- paste('{',
                    paste('"', names(x), '":', ifelse(is.na(x), 'NaN', x), sep = '', collapse = ','),
                 '}', sep = '')
  else if (length(x) > 1)
    result <- paste('[', paste(ifelse(is.na(x), 'NaN', x), collapse = ',' ), ']', sep = '')
  else
    result <- paste(ifelse(is.na(x), 'NaN', x))
  class(result) <- "json"
  result
}

json.logical <- function(x, ...) {
  if (!is.null(names(x)))
    result <- paste('{',
                 paste('"', names(x), '" : ', ifelse(x, 'true', 'false'), sep = '', collapse = ','),
                 '}', sep = '')
  else if (length(x) > 1)
    result <- paste('[', paste(ifelse(x, 'true', 'false'), collapse = ',' ), ']', sep = '')
  else
    result <- paste(ifelse(x, 'true', 'false'))
  class(result) <- "json"
  result
}

escapeStrings <- function(s) {
  s <- gsub('/', '\\\\/', s)
  s <- gsub('"', '\\\\"', s)
  s <- gsub('\\\\', '\\\\', s)
  s <- gsub('\b', '\\\\b', s)
  s <- gsub('\f', '\\\\f', s)
  s <- gsub('\n', '\\\\n', s)
  s <- gsub('\r', '\\\\r', s)
  s <- gsub('\t', '\\\\t', s)
  s
}

json.character <- function(x, ...) {
  if (!is.null(names(x)))
    result <- paste('{',
                 paste('"', names(x), '" : "', escapeStrings(x), '"', sep = '', collapse = ','),
                 '}', sep = '')
  else if (length(x) > 1)
    result <- paste('[', paste('"', escapeStrings(x), '"', collapse = ',', sep = ''), ']', sep = '')
  else
    result <- paste('"', escapeStrings(x), '"', sep = '')
  class(result) <- "json"
  result
}

json.default <- function(x, ...) {
  if (is.list(x))
    result <- json.list(x)
  else if (is.numeric(x))
    result <- json.numeric(x)
  else if (is.character(x))
    result <- json.character(x)
  else if (!is.null(names(x)))
    result <- paste('{',
                    paste('"', names(x), '" : "', escapeStrings(as.character(x)), '"', sep = '', collapse = ','),
                    '}', sep = '')
  else if (length(x) > 1)
    result <- paste('[', paste('"', escapeStrings(as.character(x)), '"', collapse = ',', sep=''), ']', sep = '')
  else
    result <- paste('"', escapeStrings(as.character(x)), '"', sep = '')
  class(result) <- "json"
  result
}

json.NULL <- function(x, ...) {
  result <- "null"
  class(result) <- "json"
  result
}

json.data.frame <- json.list
print.json <- function(x, ...) cat(x, '\n')

# This complicated regex:
# cat(gsub('((?:"[^"]*"|[^"{])*){','\\1list(', '{"a":"sdf{asdf", "b": 2}'  , perl=T) , "\n")
# from: http://regexadvice.com/forums/thread/23957.aspx
# could be used to convert json back to R
#
