add_missing = function(yr, id) { 
  ## all possible measures = years x idtrees
  dat_compl = expand.grid(year=unique(yr),idtree=unique(id))
  dat_compl$id = paste(dat_compl$idtree,dat_compl$year)
  return(dat_compl$id)
}
