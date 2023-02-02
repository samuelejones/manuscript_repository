### winner's curse function
wcc = function(BetaXG=BetaXG,seBetaXG=seBetaXG,p=5E-8){
   t = qnorm(1-p/2)
   biasA = seBetaXG*dnorm(BetaXG/seBetaXG - t)
   biasB = pnorm(BetaXG/seBetaXG  - t)
   BetaXGc = BetaXG - biasA/biasB
   # make sure winner's curse correction doesn't go below 10% of original exposure beta
   BXGc     = matrix(nrow=length(BetaXGc),ncol=2)
   BXGc[,1] = BetaXGc
   BXGc[,2] = 0.1*BetaXG
   BetaXGwc = apply(BXGc,1,max)
   return(BetaXGwc)
}


wcc2 = function(BetaXG=BetaXG,seBetaXG=seBetaXG,BetaYG=BetaYG,seBetaYG=seBetaYG,p=5E-8) {

   # Stability parameter S in below script would ideally be 1
   ### Run loop going down from 1 in intervals of 0.01 until you get convergence - to decide value of S

   Sworks <- 0

   for (stab in rev(seq(from = 0, to = 1, by = 0.01))) {
      
      break_next <- TRUE
      
      tryCatch({
      
         tryCatch({
      
            S <- stab
            #p        =  0.00000005
            t        = qnorm(1-p/2)
            L         = length(BetaXG)
            BetaXGc   = NULL
            seBetaXGc = NULL

            Wcorrect = function(a){
               gamma = a[1]
               f1 = dnorm(b,gamma,se)
               f2 = pnorm((gamma/se) - (S*t))
               l = -log(f1/f2)
            }

            for(i in 1:L){
               b            = BetaXG[i]
               se           = seBetaXG[i]
               wcresults    = optim(b, Wcorrect,method = "BFGS", hessian = TRUE)
               BetaXGc[i]   = wcresults$par
               seBetaXGc[i] = sqrt(1/wcresults$hessian)
            }

         }, error=function(e){ break_next <<- FALSE})

         #stats$BetaXG <- BetaXGc
         #stats$seBetaXG <- seBetaXGc

         # Run MR
         # Multiple outcomes of interest - use code below (will create table of results that can be used for plots)
         # Generate new variables

         tmpstats <- data.frame(BetaXG,seBetaXGc,BetaYG,seBetaYG)

         # Run heterogeneity test (if it fails, run next loop of stability check)
         p.hetero = (1-pchisq(metagen(BetaYG/BetaXG,seBetaYG/BetaXG)$Q,metagen(BetaYG/BetaXG,seBetaYG/BetaXG)$df.Q))

      }, error=function(e){ break_next <<- FALSE})

      if(break_next) { Sworks <- Sworks + 1
         # Use stab = S - 0.02 of the maximum value of S that successfully works
         if (Sworks == 3) { break } }
   }

   # where winner's curse correction either flips the sign of the exposure data
   # or attenuates to less than 10% of original, exclude variant by setting to missing
   #BetaXGc[(BetaXGc>0&BetaXG>0&BetaXGc<0.1*BetaXG)|(BetaXGc<0&BetaXG<0&BetaXGc>0.1*BetaXG)|(BetaXGc*BetaXG<0)] <- NA
   #seBetaXGc[(BetaXGc>0&BetaXG>0&BetaXGc<0.1*BetaXG)|(BetaXGc<0&BetaXG<0&BetaXGc>0.1*BetaXG)|(BetaXGc*BetaXG<0)] <- NA
   #BetaXGc[BetaXGc/BetaXG<0.1] <- NA
   #seBetaXGc[is.na(BetaXGc)] <- NA

   BetaYGc <- BetaYG
   seBetaYGc <- seBetaYG

   BetaYGc[is.na(BetaXGc)] <- NA
   seBetaYGc[is.na(BetaXGc)] <- NA

   return(data.frame(BetaXGc,seBetaXGc,BetaYGc,seBetaYGc))

}
