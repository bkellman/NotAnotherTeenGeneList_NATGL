#source('HuNet_Distance/HuNet.func.R')
library(reshape)
library(igraph)
library(RamiGO)
#library(graph)
#library(Rgraphviz)

load_HuNet <- function(threshold=.7){
	# Make a graph (g, class: graphNEL) from Human PPI db (Hippie)
	g <- Make_Protein_Graph(threshold)  

	# convert g to 'igraph' graph
	HuNet <- igraph.from.graphNEL(g, name = TRUE, weight = TRUE,unlist.attrs = TRUE)
	HuNet <- HuNet %>% set_vertex_attr("label", value = V(HuNet)$name)
	return(HuNet)
}

load_GeneMania <- function(combined=TRUE,type=c("Co-expression","Co-localization","Genetic_interactions","Pathway","Physical_interactions","Predicted","Shared_protein_domains"),alltypes=FALSE){
	if(combined){
		r <- read.csv('../data/DATA_Autism_Genomic_Varients/GeneMania/combined/COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt',sep='\t',stringsAsFactors=FALSE)
		colnames(r)[3] <- 'weight'
		g <- graph.data.frame(r,directed=FALSE)
	}else{
		ls <- system('ls ../data/DATA_Autism_Genomic_Varients/GeneMania/individual/*.txt',intern=TRUE)
		if(!any(grepl(type,ls))){stop(paste(type,'is not a valid type'))}
		rL <- list()
		for(l in ls){
			lp <- strsplit(rev(strsplit(l,split='/')[[1]])[1],split='\\.')[[1]]
			if(any(lp[1]==type) | alltypes){
				rL[[l]] <- read.csv(l,sep='\t',stringsAsFactors=FALSE)
#				rL[[l]] <- read.csv(paste('../data/DATA_Autism_Genomic_Varients/GeneMania/individual/',l,sep=''),sep='\t',stringsAsFactors=FALSE)
				rL[[l]]$type <- lp[1]
				rL[[l]]$database <- lp[2]
			}
		}
		r <- do.call(rbind,rL)
		colnames(r)[3] <- 'weight'
		g <- graph.data.frame(r,directed=FALSE)
	}
	g
}

Network_ShortestPaths_distance <- function(terms,Network,func=function(x) x,weighted=TRUE){
	terms_i = which(names(V(Network)) %in% unique(terms))
	if(weighted){
		sp <- shortest.paths(Network, v=terms_i, to=terms_i ,mode = "all",weights=NULL)
	}else{
		sp <- shortest.paths(Network, v=terms_i, to=terms_i ,mode = "all",weights=NA)
	}
	func(sp)
}

GO_jaccard<-function(BM,gene_ids,gene_type,ramigo_search=FALSE,ramigo_paths=FALSE,debug=FALSE){
	O <- list()
	for(g_i in gene_ids){
		if(debug){print(paste('get GO for',g_i))}
		o_i = BM$go_id[BM[gene_type]==g_i]
		if(nchar(paste(o_i,collapse=''))>0 ){
			if(ramigo_search){
				O[[g_i]] = GO_Ramigo_flatExpand(o_i)
			}else{
				O[[g_i]]  = o_i
			}
		}
	}
	
	sim <- matrix(0,length(gene_ids),length(gene_ids),dimnames=list(gene_ids,gene_ids))
	for(g_i in gene_ids){
		for(g_j in gene_ids){
			o_i = O[[g_i]]
			o_j = O[[g_j]]
			if(nchar(paste(o_i,collapse=''))>0 & nchar(paste(o_j,collapse=''))>0){
				if(ramigo_paths){
					sim[g_i,g_j] <- GO_Ramigo_mostRecentAncestor(o_i,o_j)
				}else{
					sim[g_i,g_j] <- length(intersect(o_i,o_j)) / length(union(o_i,o_j))
				}
			}			
		}
	}
	filter <- apply(sim,1,sum)>0 | apply(sim,2,sum)>0
	sim[filter,filter]
	sim
}

Ramigo_wrapper <- function(go,attempts=10,file='ramigo.dot'){
	dotRes <- NULL
	while(is.null(dotRes)){
		try( dotRes <- getAmigoTree(goIDs=go, picType="dot", saveResult=TRUE,filename=file) )
		attempts = attempts - 1
		if(attempts==0){stop('server cannot be reached with multiple attempts')}
	}
	#tt <- readAmigoDot(object=dotRes)
	tt <- readAmigoDot_man(file)
}

readAmigoDot_man <- function(file){
	r <- readLines('ramigo.dot')
	edges <- list()
	nodes <- list()
	for(r_i in r){
		if(grepl('->',r_i)){
			# edge
			edges[[r_i]] <- gsub('\t','',strsplit(r_i,' ')[[1]][c(1,3)])
		}else if(grepl('GO|UBERON|CARO|CL|OBA',r_i )) {
			# node
			nodes[[r_i]] <- strsplit(r_i,'>|<')[[1]][8]
		}else{
			print(paste('skip',r_i))
		}
	}
	nodes<-unlist(nodes)
	edges<-lapply(edges,function(x) unlist(nodes[as.numeric(unlist(gsub('node','',x)))]))
	eL <-do.call(rbind,edges)
	rownames(eL) <- NULL
	colnames(eL) <- c('parent','child')
	graph.data.frame(eL,directed=FALSE)
}

GO_Ramigo_flatExpand <- function(go){
	tmp<<-go
	tt <- Ramigo_wrapper(go)
	#dotRes <- getAmigoTree(goIDs=go, picType="dot", saveResult=F)
	#tt <- readAmigoDot(object=dotRes)
	#tt@annot$GO_ID
	as.character(V(tt)$name)
}


GO_Ramigo_mostRecentAncestor <- function(go1,go2,func=median){
	go <- c(go1,go2)
	tt <- Ramigo_wrapper(go)
	#dotRes <- getAmigoTree(goIDs=go, picType="dot", saveResult=F)
	#tt <- readAmigoDot(object=dotRes)
	#adjM <- tt@adjMatrix
	#rownames(adjM)=colnames(adjM)=tt@annot$GO_ID
	#graph <- graph.adjacency(adjM, mode="undirected") # shortest path in an undirected graph will give twice the distance to the "most recent common ancestor"
	graph <- tt
	sp <- shortest.paths(graph, v=go1[-(go1 %in% go2)], to=go2[-(go2 %in% go1)] ,mode = "all")
	func(sp)
}

combine <- function(dL , func,BM,gene_ids,gene_type,weight=rep(1,length(dL))){
	#initialize final distance matrix
	d <- matrix(0,length(gene_ids),length(gene_ids),dimnames=list(gene_ids,gene_ids))
	
	# convert to gene_type names
	for(i in 1:length(dL)){
		d_i = dL[[i]]
		d_i_id <- names(which.max(apply(BM,2,function(x) sum(rownames(d_i) %in% x))))
		if(!all(rownames(d_i)==colnames(d_i))){stop('distance matrix is not symetrical')}
		rownames(d_i)=colnames(d_i)=convert(rownames(d_i),d_i_id,gene_type,BM)
		dL[[i]] = d_i
	}
	# apply weights
	for(i in 1:length(dL)){
		dL[[i]] = dL[[i]] * weight[i]
	}
	# combine distance matrixes to d
	for(i in gene_ids){
		for(j in gene_ids){
			d[i,j] = func(na.omit(as.vector(unlist(lapply(dL,function(x) get_dij(x,i,j))))))
		}
	}
	return(d)
}

get_dij <- function(d,i,j){
	if(i %in% rownames(d) & j %in% colnames(d)){
		d[i,j]
	}else{
		median(d)
		print(paste(i,'or',j,'not found in d'))
	}
}

convert <- function(a,aname,bname,BM){
	unlist(sapply(a,function(a_i) unique( BM[[bname]][which(BM[[aname]]==a_i)])))
}





