
# author : Julien Tap
# date : 2016 April 18th
# licence : GPL-3


# rat_mixer is a function to mix rat individual according 3 constrains : weights, microbiota and filiations.
# 


rat_mixer =  function(weight, microbiota, filiation, cage_provider=NULL,  groups_nb=7, nb_filiation_max = 2, nb_cage_max=5,  microbiota_quantile=TRUE, plot=FALSE ){

if(!is.numeric(weight)) stop("weight have to be numeric")
if(!is.numeric(microbiota)) stop("microbiota have to be numeric")
if(!is.factor(filiation)) stop("filiation have to be a factor")

if(!is.numeric(groups_nb)) stop("groups_nb have to be numeric")
l = c(length(weight), length(microbiota), length(filiation))
if(length(unique(l)) != 1) stop("weight, microbiota, filiation do not have the same length")

if(!is.null(cage_provider)){

	if(!is.factor(cage_provider)) stop("cage_provider have to be a factor")

}


l=l[1] # nb of individuals

create_folds = function(l, nfolds) {

	perm = sample(1:l, l) / l
	id   = rep(0, l)

	for (f in nfolds:1) { id[perm <= f/nfolds] = f }

	return(id)

}

balanced_all = FALSE
i = 1

while(!balanced_all){
	
	groups = create_folds(l=l,nfolds=groups_nb)
	
	if(microbiota_quantile) {
	
		groups_nb2 = round(mean(table(groups)))
		groups_new = create_folds(l=l,nfolds=groups_nb2)
		groups = sort(groups_new)[order(order(microbiota))]
		groups_new=rep(NA, length(groups))
	
		for(j in 1:groups_nb2){
			
			groups_new[which(groups==j)] = sample(1:length(which(groups==j)))
			#groups_new[which(groups==j)] = 1:length(which(groups==j))
		}
		
		groups = groups_new
	}
	
	
	p_weight     = na.omit(suppressWarnings(as.vector(pairwise.wilcox.test(weight,groups, p.adjust="none")$p.value)))
	p_microbiota = na.omit(suppressWarnings(as.vector(pairwise.wilcox.test(microbiota,groups, p.adjust="none")$p.value)))

	balanced_weight     = length(which(c(p_weight > 0.05) == FALSE)) == 0
	balanced_microbiota = length(which(c(p_microbiota > 0.05) == FALSE)) == 0
	balanced_filiation  = length(which(as.vector(apply(table(groups, filiation),1, function(x) x <= nb_filiation_max )) == FALSE)) == 0

	if(!is.null(cage_provider)){
	
	balanced_cage_provider  = length(which(as.vector(apply(table(groups, cage_provider),1, function(x) x <= nb_cage_max )) == FALSE)) == 0
	
	}
	
	if(is.null(cage_provider)) {
	balanced_all = length(which(c(balanced_weight, balanced_microbiota, balanced_filiation) == FALSE)) == 0
	
	} else { balanced_all = length(which(c(balanced_weight, balanced_microbiota, balanced_filiation, balanced_cage_provider) == FALSE)) == 0 }
	

		i=i+1
		
		if (i>99) stop("your constrains are too high, increase nb_filliation_max for example")

	}
	
cat(i, " iterations were needed to balance the dataset\n")

if(plot){
	par(mfrow=c(2,2))
	boxplot(weight~groups, ylab="weight", xlab="groups")
	boxplot(microbiota~groups, ylab="microbiota input", xlab="groups")
	barplot(table(filiation,groups), col=rainbow(length(levels(filiation))), xlab="groups", main="filiation")
	
	if(!is.null(cage_provider)) {
	barplot(table(cage_provider,groups), col=rainbow(length(levels(cage_provider))), xlab="groups", main="cage provider")
	}
	
}

return(groups)



}
