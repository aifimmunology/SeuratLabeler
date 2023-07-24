#' Write a stderr message with a leading date/time stamp
#'
#' @param x a character object with the message to display
#'
#' @return no return
#' @export
stm <- function(x) {
  assertthat::assert_that(class(x) == "character")
  assertthat::assert_that(length(x) == 1)

  write(paste0("[",Sys.time(),"] ",x), stderr())
}

#' Suppress selected messages and/or warnings from function output
#' 
#' Returns the result slot of the quietly call. Prints all output in the outputs slot.
#' Outputs messages for all messages in the messages slot that do not match expected
#' values or patterns. Outputs warnings for all items in warnings slot that do not
#' match expected values or patterns.
#'
#' @param A function call result generated using \code{purrr::quietly()}
#' @param values_suppress Named list of character vectors. Names should be be 'messages' and 'warnings'.
#' Vectors are values of expected messages and warnings that should be suppressed (exact matches). For 
#' example, \code{list(messages=c('Expected message'), warnings = c("Expected warning"))}
#' @param patterns_suppress Named list of character vectors. Names should be be 'messages' and 'warnings'.
#' Vectors are patterns of expected messages and warnings that should be suppressed (pattern matched)
#' @returns The $result value of the input quietly call result.
#' @examples
#' testfun <- function(){
#'     print('testprint1')
#'     print('asdfsdaf')
#'     print('sdfasdfas')
#'     message('expected_message_1')
#'     message('exp_message_2')
#'     message('unexpected message1')
#'     warning('exp_warn_2')
#'     warning('expected warn')
#'     warning('unexpected warn')
#'     return(NULL)
#' }
#' vlist <- list(messages=c('expected_message_1'), warnings = c("exp_warn_2"))
#' plist <- list(messages=c('test'), warnings = c("^expected"))
#' # # See result of filtered quietly call
#' # filter_quietly(purrr::quietly(testfun)(), values_suppress=vlist, plist)
#' # # Compare to result of original function call
#' # testfun()
filter_quietly <- function(
    quietly_res, 
    values_suppress = list(messages = character(), 
                           warnings = character()),
    patterns_suppress = list(messages = character(), 
                             warnings = character())
){
    default_list <- list(
        messages = character(), 
        warnings = character()
    )

    # Validate values_suppress
    missing_val <- setdiff(names(default_list), names(values_suppress))
    extra_val <- setdiff(names(values_suppress), names(default_list))
    if (length(missing_val) > 0){
        values_suppress <- c(values_suppress, default_list[missing_val])
    }
    if (length(extra_val) > 0){
        warning(sprintf("Unexpected list item names in values_suppress: %s. Expected labels are messages, output, or warnings",
                       paste(extra_val, collapse = ',')))
    }
    
    # Validate patterns_suppress
    missing_pat <- setdiff(names(default_list), names(patterns_suppress))
    extra_pat <- setdiff(names(patterns_suppress), names(default_list))
    if (length(missing_pat) > 0){
        patterns_suppress <- c(patterns_suppress, default_list[missing_pat])
    }
    if (length(extra_pat) > 0){
        warning(sprintf("Unexpected list item names in patterns_suppress: %s. Expected labels are messages, output, or warnings",
                       paste(extra_pat, collapse = ',')))
    }
    
    .match_listwise = function(query_list, pattern_list){
        matchvals <- unlist(sapply(pattern_list, function(x){grep(x, query_list, value = T)}))
        if(length(matchvals)>0){
            return(matchvals)
        } 
        
    }
    .filter_values <- function(fieldname, patterns_suppress, values_suppress){
        ignore_match <- .match_listwise(quietly_res[[fieldname]], patterns_suppress[[fieldname]])
        ignore_value_i <- which(gsub("\n","", quietly_res[[fieldname]]) %in% values_suppress[[fieldname]])
        ignore_value <- quietly_res[[fieldname]][ignore_value_i]
        ignore_all <- unique(c(ignore_match, ignore_value))
        output_vals <- setdiff(quietly_res[[fieldname]], ignore_all)
        return(output_vals)
    }

    # Print all output
    if(length(quietly_res$output) > 0){
        for (outmsg in quietly_res$output) {
            cat(outmsg)
        }
    }
    
    # Handle messages
    message_vals <- .filter_values('messages', patterns_suppress, values_suppress)
    if(length(message_vals) > 0){
        for (msg in message_vals) {
            message(msg)
        }
    }
    
    # Handle warnings
    warn_vals <- .filter_values('warnings', patterns_suppress, values_suppress)
    if(length(warn_vals) > 0){
        for (warnval in warn_vals) {
            warning(warnval, call.= FALSE)
        }
    }    
    
    return(quietly_res$result)
    
}