carefulFit=function(model,traitExtension,manifestTraitvarExtension, weighting = 100){

    #     if(transformedParams==TRUE) model$negDRIFTlog$free[row(DRIFT$free)!=col(DRIFT$free)] <- FALSE #fix off diagonal DRIFT params
    #     if(transformedParams==FALSE) model$DRIFT$free[row(DRIFT$free)!=col(DRIFT$free)] <- FALSE
    
    if(traitExtension==TRUE) penalties <- OpenMx::mxAlgebra(name='penalties',
      sum(T0VAR*T0VAR) + sum(1/((diag2vec(T0VAR) * diag2vec(T0VAR)))) +
        sum(DRIFT*DRIFT) + sum(1/(diag2vec(DRIFT) * diag2vec(DRIFT))) + sum(diag2vec(DRIFT)) +
        sum(DIFFUSION*DIFFUSION) + sum(1/((diag2vec(DIFFUSION) * diag2vec(DIFFUSION)))) +
        sum(MANIFESTVAR*MANIFESTVAR) + sum(1/((diag2vec(MANIFESTVAR) * diag2vec(MANIFESTVAR)))) +
        sum(TRAITVAR * TRAITVAR) + sum(1/((diag2vec(TRAITVAR) * diag2vec(TRAITVAR)))) +
        sum(T0TRAITEFFECT * T0TRAITEFFECT) - sum((diag2vec(T0TRAITEFFECT) * diag2vec(T0TRAITEFFECT)))
    )

    if(manifestTraitvarExtension==TRUE) penalties <- OpenMx::mxAlgebra(name='penalties',
      sum(T0VAR*T0VAR) + sum(1/((diag2vec(T0VAR) * diag2vec(T0VAR)))) +
        sum(DRIFT*DRIFT) + sum(1/(diag2vec(DRIFT) * diag2vec(DRIFT))) + sum(diag2vec(DRIFT)) +
        sum(DIFFUSION*DIFFUSION) + sum(1/((diag2vec(DIFFUSION) * diag2vec(DIFFUSION)))) +
        sum(MANIFESTVAR*MANIFESTVAR) + sum(1/((diag2vec(MANIFESTVAR) * diag2vec(MANIFESTVAR)))) +
        sum(MANIFESTTRAITVAR * MANIFESTTRAITVAR) + sum(1/((diag2vec(MANIFESTTRAITVAR) * diag2vec(MANIFESTTRAITVAR))))
    )

    if(traitExtension==FALSE & manifestTraitvarExtension==FALSE)  penalties <- OpenMx::mxAlgebra(name='penalties',
      sum(T0VAR*T0VAR) + sum(1/((diag2vec(T0VAR) * diag2vec(T0VAR)))) +
        sum(DRIFT*DRIFT) + sum(1/(diag2vec(DRIFT) * diag2vec(DRIFT))) + sum(diag2vec(DRIFT)) +
        sum(DIFFUSION*DIFFUSION) + sum(1/((diag2vec(DIFFUSION) * diag2vec(DIFFUSION)))) +
        sum(MANIFESTVAR*MANIFESTVAR) + sum(1/((diag2vec(MANIFESTVAR) * diag2vec(MANIFESTVAR))))
    )
  
  # penalties <- OpenMx::mxAlgebra(name='penalties',
  #       sum(DRIFT*DRIFT) + sum(1/(diag2vec(DRIFT) * diag2vec(DRIFT)))   
  #   )
    
    penaltyLL <- OpenMx::mxAlgebra(sum(ctsem.fitfunction)+ctsem.penalties*FIMLpenaltyweight, name='penaltyLL')
    
    modelwithpenalties <- OpenMx::mxModel(model, 
      penalties, 
      mxFitFunctionML(vector=FALSE)
    )
    
    model<-OpenMx::mxModel('ctsemCarefulFit', 
      modelwithpenalties, penaltyLL,
      #             mxMatrix(type='Full', name='FIMLpenaltyweight', nrow=1, ncol=1, values=FIMLpenaltyweight, free=F), 
      mxMatrix(name='FIMLpenaltyweight', values=weighting,free=FALSE,nrow=1,ncol=1,type='Full' ), 
      mxFitFunctionAlgebra('penaltyLL')
    )
    
    model<-try(suppressWarnings(OpenMx::mxRun(model)))
  }
