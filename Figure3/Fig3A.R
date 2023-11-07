P=SpatialFeaturePlot(object = sample2_norm, 
                   features = c("MS4A1"), 
                   alpha = c(0.5, 1),pt.size.factor = 1)
P+ggplot2::scale_fill_gradient2(low="green", mid="yellow", high="red",midpoint = 0.5,na.value = "red",space="lab",limits = c(0,1),breaks =c(0,0.5,1))
