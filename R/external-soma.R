soma = function (costFunction, bounds, options = list(), init = NULL, ...) {
  if (!inherits(bounds, "soma.bounds")) 
    bounds = do.call(soma::bounds, bounds)
  if (!inherits(options, "soma.options")) 
    options = do.call(soma::all2one, options)
  nParams = length(bounds$min)
  nParamsTotal = nParams * options$populationSize
  if (!is.null(options$pathLength) && !is.null(options$stepLength)) 
    steps = seq(0, options$pathLength, options$stepLength)
  nSteps = options$nSteps %||% length(steps)
  perturbationChance = options$perturbationChance %||% 1
  leaderPoolIndices = migrantPoolIndices = NULL
  if (options$strategy == "pareto") {
    leaderPoolIndices = seq_len(ceiling(options$populationSize/25))
    migrantPoolIndices = seq_len(ceiling(options$populationSize/6.25)) + 
      round(options$populationSize/5)
  }
  if (!is.null(init)) 
    population = init
  else {
    population = matrix(runif(nParamsTotal), nrow = nParams, 
                         ncol = options$populationSize)
    population = population * (bounds$max - bounds$min) + 
      bounds$min
  }
  evaluations = 0
  evaluationCount = 0
  costFunctionWrapper = function(x) {
    evaluationCount <= evaluationCount + 1
    costFunction(x, ...)
  }
  costFunctionValues = apply(population, 2, costFunctionWrapper)
  migrationCount = 0
  leaderCostHistory = min(costFunctionValues)
  # report(OL$Info, "Starting SOMA optimisation")
  repeat {
    progress = migrationCount/options$nMigrations
    if (options$strategy == "all2one") {
      leader = which.min(costFunctionValues)
      migrants = seq_len(options$populationSize)
    }
    else if (options$strategy == "t3a") {
      perturbationChance = 0.05 + 0.9 * progress
      steps = (seq_len(nSteps) - 1) * (0.15 - 0.08 * 
                                          progress)
      leaderPool = sample(options$populationSize, options$leaderPoolSize)
      leader = leaderPool[which.min(costFunctionValues[leaderPool])]
      migrantPool = sample(options$populationSize, options$migrantPoolSize)
      migrants = migrantPool[order(costFunctionValues[migrantPool])[seq_len(options$nMigrants)]]
    }
    else if (options$strategy == "pareto") {
      perturbationChance = 0.5 + 0.45 * cos(options$perturbationFrequency * 
                                               pi * progress + pi)
      steps = (seq_len(nSteps) - 1) * (0.35 + 0.15 * 
                                          cos(options$stepFrequency * pi * progress))
      order = order(costFunctionValues)
      leader = sample(order[leaderPoolIndices], 1)
      migrants = sample(order[migrantPoolIndices], 1)
    }
    leaderValue = costFunctionValues[leader]
    separationOfExtremes = max(costFunctionValues) - min(costFunctionValues)
    sumOfExtremes = max(costFunctionValues) + min(costFunctionValues)
    if (migrationCount == options$nMigrations) {
      # report(OL$Info, "Migration limit (#{options$nMigrations}) reached - stopping")
      break
    }
    if (separationOfExtremes < options$minAbsoluteSep) {
      # report(OL$Info, "Absolute cost separation (#{separationOfExtremes}) is below threshold (#{options$minAbsoluteSep}) - stopping", 
             # signif = 3)
      break
    }
    if (isTRUE(abs(separationOfExtremes/sumOfExtremes) < 
               options$minRelativeSep)) {
      # report(OL$Info, "Relative cost separation (#{separationOfExtremes/sumOfExtremes}) is below threshold (#{options$minRelativeSep}) - stopping", 
             # signif = 3)
      break
    }
    toPerturb = array(FALSE, dim = dim(population))
    toPerturb[, migrants] = runif(nParams * length(migrants)) < 
      perturbationChance
    if (options$strategy == "pareto") 
      toPerturb[, migrants] = ifelse(toPerturb[, migrants], 
                                      1, progress)
    toMigrate = colSums(toPerturb) > 0
    toMigrate[leader] = FALSE
    nMigrating = sum(toMigrate)
    if (nMigrating == 0) {
      # report(OL$Verbose, "No parameters to perturb - skipping migration")
      next
    }
    # else report(OL$Verbose, "Migration ##{migrationCount+1}: #{nMigrating} individuals moving towards leader ##{leader} with cost #{leaderValue}", 
    #             signif = 3)
    directionsFromLeader = population[, toMigrate] - population[, 
                                                                 leader]
    populationSteps = array(rep(population[, toMigrate], 
                                 nSteps), dim = c(nParams, nMigrating, nSteps))
    populationSteps = populationSteps - rep(steps, each = nParams * 
                                               nMigrating) * rep(directionsFromLeader * toPerturb[, 
                                                                                                  toMigrate], nSteps)
    outOfBounds = which(populationSteps < bounds$min | 
                           populationSteps > bounds$max)
    randomSteps = array(runif(nParamsTotal * nSteps), dim = c(nParams, 
                                                               options$populationSize, nSteps))
    randomSteps = randomSteps * (bounds$max - bounds$min) + 
      bounds$min
    populationSteps[outOfBounds] = randomSteps[outOfBounds]
    allCostFunctionValues = apply(populationSteps, 2:3, 
                                   costFunctionWrapper)
    individualBestLocs = apply(allCostFunctionValues, 1, 
                                which.min)
    indexingMatrix = cbind(seq_len(nMigrating), individualBestLocs)
    subpopulation = t(apply(populationSteps, 1, "[", indexingMatrix))
    population[, toMigrate] = subpopulation
    bestCostFunctionValues = allCostFunctionValues[indexingMatrix]
    costFunctionValues[toMigrate] = bestCostFunctionValues
    migrationCount = migrationCount + 1
    evaluations = c(evaluations, evaluationCount)
    leaderCostHistory = c(leaderCostHistory, min(costFunctionValues))
  }
  leader = which.min(costFunctionValues)
  # report(OL$Info, "Final leader is ##{leader}, with cost #{costFunctionValues[leader]}", 
         # signif = 3)
  returnValue = list(leader = leader, population = population, 
                      cost = costFunctionValues, history = leaderCostHistory, 
                      migrations = migrationCount, evaluations = evaluations)
  class(returnValue) = "soma"
  return(returnValue)
}

"%||%" <- function (X, Y) { if (is.null(X) || length(X)==0) Y else X }

