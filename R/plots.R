#' Plot recursive segmentation tree
#'
#' Plot the tree resulting from a recursive segmentation procedure
#'
#' @param tree data frame, tree resulting from the call of Segmentation_Recursive
#'
#' @return a ggplot
#'
#' @examples
#' res=Segmentation_Recursive(obs=RhoneRiverAMAX$H)
#' plot(res$data$obs,col=res$data$period)
#' PlotTree(res$tree)
#' @export
#' @import ggplot2
#' @importFrom stats rnorm setNames
#' @importFrom scales viridis_pal
PlotTree <- function(tree){
  DF=tree
  n=NROW(DF)
  if(n<=1){
    message('No shift detected, no tree structure to plot')
    return(NULL)
  }
  # Add columns for x/y's of nodes and arrows
  DF$y=-1*DF$level
  DF$x=0*DF$y
  DF$xstart=DF$xend=DF$ystart=DF$yend=NA*DF$y
  DF$isTerminal=FALSE

  # Assign x coordinates for all nodes in a same level
  x_ref_plot=0 # central axis fixed at 0 for the first level
  for(i in 2:max(DF$level)){# first level skipped, always x = 0
    indx_table=which(DF$level==i)
    # Get the number of nodes in this level and assign coordinates with middle point at x=0
    half_n <- (length(indx_table) - 1) / 2
    coord_local_node <- seq(-half_n, half_n)

    id_parents=unique(DF$parent[indx_table])

    # if same parent -> reinitialize x offset
    if(length(id_parents)==1){
      x_ref_plot = DF$x[id_parents]
    }
    # Allocate coordinates of nodes at this level adding x_ref_plot estimated of the last level (offset) to move the reference axis
    DF$x[indx_table] = coord_local_node+x_ref_plot

    # Give a weight to each node in function of the number of children. This helps to offset the reference axis to plot new nodes
    # There are more children on the left than on the right, so moving to the left makes more sense for the plot of the next level.
    weight_x_plot =  ifelse(coord_local_node+x_ref_plot < 0, -1, ifelse(coord_local_node+x_ref_plot > 0, 1, 0)) *  DF$nS[indx_table]/max(DF$nS[indx_table])

    x_ref_plot=x_ref_plot+
      (sum(weight_x_plot)) # move the xref_plot to plot nodes in this level

  }
  # Assign segments to link nodes
  if(n>1){
    for(i in 2:n){
      # arrows
      DF$xstart[i]=DF$x[DF$parent[i]]
      DF$xend[i]=DF$x[i]
      DF$ystart[i]=DF$y[DF$parent[i]]-0.1
      DF$yend[i]=DF$y[i]+0.1
      # is node terminal?
      DF$isTerminal[i]=sum(DF$parent==DF$indx[i])==0
    }
  } else {
    DF$isTerminal[1]=TRUE
  }

  DF$fontface='plain';DF$fontface[DF$isTerminal]='bold'
  # Define the colors and remove the level 0
  tree_levels <- unique(DF$parent)[-1]
  tree_colors <-scales::viridis_pal(option='D')(length(tree_levels))

  # Create a named vector of colors for the legend
  color_legend <- stats::setNames(tree_colors, as.character(tree_levels))
  # Get the first node
  DF_initial_node=DF[1,]
  # Get tree structure without first node
  DF_tree_plot=DF[-1,]

  # plot
  g=ggplot(DF_tree_plot)+
    # plot initial node
    geom_point(data=DF_initial_node,
               aes(x,y),
               size=11, col='gray', shape=16, alpha=0.5)+
    # Text for the initial node
    geom_text(data=DF_initial_node,
              aes(x,y,label=indx),
              size=4)+
    # Arrows to link the nodes
    geom_segment(aes(x=xstart,y=ystart,xend=xend,yend=yend),
                 arrow=arrow(type='closed',length=unit(0.015, "npc")))+
    # Plot nodes of the tree structure
    geom_point(aes(x,y,col=factor(parent),shape=factor(isTerminal)),
               size=8,alpha=0.8)+
    geom_text(aes(x,y,label=indx,fontface=fontface),
              size=4)+
    scale_color_manual(name='Parent node',
                       values=color_legend,
                       breaks=tree_levels,
                       labels = as.character(tree_levels))+
    scale_shape_manual('Terminal node',values=c(15,17))+
    theme_void()+
    theme(plot.margin = margin(20, 20, 20, 20),
          legend.key.size = unit(1, "cm"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=10))

  return(g)
}
