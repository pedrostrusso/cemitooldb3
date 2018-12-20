library(CEMiTool)
library(data.table)
library(RMySQL)
library(pool)
library(plyr)

user_data <- get_user_data("prusso", "cemitooldb3", "127.0.0.1")
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- CEMiTool::read_gmt(gmt_fname)
## Get example interactions file
int_df <- read.delim(system.file("extdata", "interactions.tsv", package = "CEMiTool"))

p <- list()
p$interactions <- int_df
p$gmt <- gmt_in

studies_to_db <- function(user_data,
                          min_sample_num=18, max_sample_num=200,
                          input_dir="/home/prusso/Documents/CSBL/cemitooldb3/data",
                          results_dir="/home/prusso/Documents/CSBL/cemitooldb/intermediate/CEMiTool",
                          verbose=TRUE, force=FALSE, ...){
    p <- list(...)
    p$verbose <- verbose
    results_dir <- file.path(results_dir, gpl_id)
    
    if(dir.exists(results_dir)){
        if(!force){
            stop("Stopping analysis: ", results_dir, 
                 " already exists! Use force=TRUE to overwrite")
        }
    }else{
        dir.create(results_dir, recursive=TRUE)
    }
    
    if(verbose) message("Creating database link...")
    pool <- pool::dbPool(drv=RMySQL::MySQL(),
                         dbname=user_data$dbname,
                         host=user_data$host,
                         user=user_data$user,
                         password=user_data$password)
    
    
    for(type in dir(input_dir)){
        message("Inputting data for ", type, ":")
        path <- file.path(input_dir, type)
        for(system_organ in dir(path)){
            message("...", system_organ)
            sys_org_path <- file.path(path, system_organ)
            for(cell_tissue in dir(sys_org_path)){
                if(file_test("-d", file.path(sys_org_path, cell_tissue))){
                    message("......", cell_tissue)
                    cell_type_path <- file.path(sys_org_path, cell_tissue)
                    for(cell_file in dir(cell_type_path)){
                        cell <- gsub("^human_([^.]+)_expression.csv", "\\1", cell_file)
                        print(cell)
                    }
                }else{
                    cell_file <- cell_tissue
                    cell <- gsub("^human_([^.]+)_expression.csv", "\\1", cell_file)
                    print(cell)
                }
            }
        }
    }
    
    for(gds_id in gds_studies){
        
        if(verbose) message("Study ", which(gds_studies == gds_id), " of ", length(gds_studies))
        if(verbose) message("Getting data from GEO for study ", gds_id)
        possible_error <- tryCatch(gds <- get_gds(gds_id),
                                   error=function(e) e)
        if(inherits(possible_error, "error")){
            message("Error in GEOquery")
            next
        }else{
            gds <- possible_error
        }
        if(verbose) message("Getting and processing expression and annotation data")
        expr <- get_expr_from_gds(gds)
        
        sample_num <- ncol(expr) - 1
        if(sample_num > max_sample_num | sample_num < min_sample_num){
            message(paste("Sample size for study", gds_id,
                          "is out of the range defined by max_sample_num and min_sample_num"))
            next
        }
        
        expr <- process_and_collapse(expr)
        annot <- make_gds_annot(gds)
        #names_list[[gds_id]] <- annot
        
        run_date <- lubridate::now()
        if(verbose) message("Running CEMiTool")
        cem <- run_cemitool(gds_id, expr, annot, results_dir, p)
        
        if(!is.null(cem)){
            if(verbose) message("Populating MySQL database")
            populate_sql(pool=pool, gds_id=gds_id, user_data=user_data, run_date=run_date, cem=cem)
        }
    }
    poolClose(pool)
}

for(gpl_id in gpl_ids2){
    studies_to_db(user_data=user_data,
                  gpl_id=gpl_id,
                  min_sample_num=18, 
                  max_sample_num=200,
                  results_dir="/home/prusso/Documents/CSBL/new_database/intermediate/test",
                  verbose=TRUE, force=FALSE)
}


# mapply(studies_to_db, gpl_id=gpl_ids2, 
#        MoreArgs = list(user_data=user_data,
#        min_sample_num=18, max_sample_num=200, 
#        results_dir="/home/prusso/Documents/CSBL/new_database/intermediate/test",
#        verbose=TRUE, force=FALSE))

