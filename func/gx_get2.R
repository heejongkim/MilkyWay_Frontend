########### gx_get2 : using wget to download the dataset
# Fixed by Hee Jong Kim
gx_get2 <- function(file_id,create=FALSE,force=FALSE){
  check_url_and_key()
  file_path = file.path(gx_get_import_directory(create=create), file_id)
  if( !force && file.exists(file_path)){
    message("You already downloaded this file, use force=TRUE to overwrite")
    return(file_path)
  }
  hist_datasets <- gx_list_history_datasets()
  encoded_dataset_id <- hist_datasets[hist_datasets$hid==file_id,'id']
  dataset_details <- gx_show_dataset(encoded_dataset_id)

  if( dataset_details$state == 'ok' ){
    url <- paste0(
      pkg.env$GX_URL,'api/histories/',pkg.env$GX_HISTORY_ID,
      '/contents/',encoded_dataset_id,'/display',
      '?to_ext=',dataset_details$extension,
      '&key=',pkg.env$GX_API_KEY)
    download.file(url,file_path,quiet=TRUE, method='wget')
  }
  return(file_path)
}
########### /gx_get2 : using wget to download the dataset
