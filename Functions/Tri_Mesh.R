Tri_Mesh <- function(m, b = NA) {
  options(warn = -1)
  #### Get points
  p_mat <- matrix(0, m, m)
  p_mat[upper.tri(p_mat, diag = T)] <- 1
  
  if(!is.na(b)){
    
    for (i in 1:m) {
      for (j in 1:m) {
        if (j - i >= b) {
          p_mat[i, j] <- 0
        }
      }
    }
    
  }
  
  
  pts <- melt(p_mat)
  
  pts <- dplyr::filter(pts, value == 1)
  pts$n <- 1:nrow(pts)
  
  
  pts <- as.matrix(pts[,1:2])
  
  
  ### Get triangles
  tri_mat <- matrix(0, 1, 3)
  
  
  for (k in 1:m-1) {
    
    
    n = which(pts[,1] == k)
    p = which(pts[,1] == k+1)
    
    
    
    for (j in 1:length(p)) {
      
      tri_mat <- rbind(tri_mat, c(n[j], n[j+1], p[j]))
      tri_mat <- rbind(tri_mat, c(p[j-1], n[j], p[j]))
      
    }
    
  }
  
  # filter out incorrect triangles
  
  rem <- c(which(tri_mat[,1] == tri_mat[,3]))
  tri_mat <- tri_mat[-rem,]
  tri_mat <- na.omit(tri_mat)
  
  ### Get edges
  
  
  edg_mat <- matrix(0, 1, 2)
  
  v <- which(pts[,1] == 1)
  h <- which(pts[,2] == m)
  d <- c()
  m_rev <- rev(c(1:(m)))
  for (i in m_rev) {
    d <- append(d, values = which(pts[,1] == i & pts[,2] == i))
  }
  
  
  for (i in 1:m) {
    ve <- c(v[i], v[i+1])
    edg_mat <- rbind(edg_mat, ve)
    
  }
  
  edg_mat <- edg_mat[-1,] %>% na.omit()
  
  for (i in 1:m) {
    he <- c(h[i], h[i+1])
    edg_mat <- rbind(edg_mat, he)
    
  }
  
  edg_mat <- na.omit(edg_mat)
  
  for (i in 1:m) {
    de <- c(d[i], d[i+1])
    edg_mat <- rbind(edg_mat, de)
    
  }
  
  edg_mat <- na.omit(edg_mat)
  pts <- pts - 1
  
  return(list("Nodes" = pts, "Triangles" = tri_mat, "Edges" = edg_mat))
  
}

