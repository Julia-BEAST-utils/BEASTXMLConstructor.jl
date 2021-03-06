module BeastNames

const TAXA = "taxa"
const TAXON = "taxon"
const ID = "id"
const NAME = "name"
const MISSING_VAL = "NA"
const NEWICK = "newick"
const MCMC = "mcmc"
const LOG = "log"
const FILENAME = "fileName"
const CHAINLENGTH = "chainLength"
const LOGEVERY = "logEvery"
const DATE = "date"
const OPERATORS = "operators"
const ATTR = "attr"
const MULTIVARIATE_DIFFUSION_MODEL = "multivariateDiffusionModel"
const PRECISION_MATRIX = "precisionMatrix"
const MATRIX_PARAMETER = "matrixParameter"
const PARAMETER = "parameter"
const VALUE = "value"
const SAMPLING_PRECISION = "samplingPrecision"
const REPEATED_MEASURES = "repeatedMeasuresModel"
const WEIGHT = "weight"
const BEAST = "beast"

const TRAIT_VALIDATION = "traitValidation"
const LOG_SUM = "logSum"
const TRAIT_VALIDATION_PROVIDER = "traitValidationProvider"
const CROSS_VALIDATION = "crossValidation"

const TRAIT_NAME = "traitName"
const INFERRED_TRAIT = "inferredTrait"
const USE_TREE_TRAITS = "useTreeTraits"
const TRAIT_DATA_LIKELIHOOD = "traitDataLikelihood"
const IDREF = "idref"
const TRAIT_PARAMETER = "traitParameter"
const FILELOG = "fileLog"
const MULT_DIST_LIKELIHOOD = "multivariateDistributionLikelihood"
const DISTRIBUTION = "distribution"
const MVM_DISTRIBUTION = "multivariateNormalDistributionModel"
const MEAN_PARAMETER = "meanParameter"
const PRECISION_PARAMETER = "precisionParameter"
const DATA = "data"
const LEAF_TRAIT_PARAMETER = "leafTraitParameter"
const TREE_MODEL = "treeModel"
const LIKELIHOOD = "likelihood"
const POSTERIOR = "posterior"
const TRUE = "true"
const FALSE = "false"
const FAST_MATRIX_PARAMETER = "fastMatrixParameter"
const ROWS = "rows"
const COLUMNS = "columns"

const DEFAULT_TRAIT_NAME = "traits"
const DEFAULT_TREE_NAME = "startingTree"

const FIXED_TREE = "fixedTree"
const TREE = "tree"
const ROOT_HEIGHT = "rootHeight"
const INTERNAL_NODES = "internalNodes"
const ROOT_NODE = "rootNode"
const NODE_HEIGHTS = "nodeHeights"
const INTERNAL_NODE_HEIGHTS = "internalNodeHeights"
const ALL_INTENAL_NODE_HEIGHTS = "allInternalNodeHeights"
const NODE_TRAITS = "nodeTraits"
const LEAF_NODES = "leafNodes"
const TRAIT_DIM = "traitDimension"
const AS_MATRIX = "asMatrix"
const LEAF_TRAITS = "leafTraits"

const MULTIVARIATE_WISHART_PRIOR = "multivariateWishartPrior"
const SCALE_MATRIX = "scaleMatrix"
const DF = "df"
const DEFAULT_MBD_PRIOR = "precisionPrior"

const DIFFUSION_ID = "diffusionModel"
const DIAGONAL_MATRIX = "DiagonalMatrix"
const DIFF_PREC_ID = "diffusion.precision"
const DIMENSION = "dimension"
const LOWER = "lower"
const UPPER = "upper"

const DEFAULT_RM_NAME = "repeatedMeasures"
const DEFAULT_RM_PREC_NAME = "residualPrecision"
const DEFAULT_RM_PREC_PRIOR_NAME = "residualPrior"

const CACHE_BRANCHES = "cacheBranches"
const ALLOW_IDENTICAL = "allowIdentical"
const USE_TREE_LENGTH = "useTreeLength"
const SCALE_BY_TIME = "scaleByTime"
const REPORT_AS_MULTIVARIATE = "reportAsMultivariate"
const INTEGRATE_INTERNAL_TRAITS = "integrateInternalTraits"
const FACTOR_TRAIT_NAME = "factors"
const CONJUGATE_ROOT_PRIOR = "conjugateRootPrior"
const PSS = "priorSampleSize"
const TRAIT_DATA_LIKELIHOOD_ID = "traitLikelihood"
const PSS_VAL = 0.001
const ALLOW_SINGULAR = "allowSingular"
const STANDARDIZE = "standardize"

const WISHART_STATISTICS = "wishartStatistics"
const REPEATED_MEASURES_WISHART_STATISTICS = "RepeatedMeasuresWishartStatistics"
const PRECISION_GIBBS_OPERATOR = "precisionGibbsOperator"
const COMPOUND_PRECISION_OPERATOR = "compoundPrecisionOperator"
const DIFFUSION_OPERATOR = "diffusionOperator"
const RESIDUAL_OPERATOR = "residualOperator"

const PRIOR = "prior"
const PRIOR_ID = "prior"
const LIKELIHOOD_ID = "likelihood"
const SCREEN_LOGEVERY = 1000
const FILE_LOGEVERY = 1000
const SCREEN_LOG_ID = "screenLog"
const FILE_LOG_ID = "fileLog"
const OVERWRITE = "overwrite"
const AUTO_OPTIMIZE = "autoOptimize"
const DEFAULT_FILENAME = "fileLog"

const COLUMN = "column"
const LABEL = "label"
const DP = "dp"
const WIDTH = "width"
const DP_DEFAULT = 4
const WIDTH_DEFAULT = 12

const TIMER = "timer"
const REPORT = "report"
const PROPERTY = "property"

const MATRIX_INVERSE = "matrixInverse"

const DIRECTION = "direction"
const FORWARDS = "forwards"

const USING_HEIGHTS = "usingHeights"
const USING_DATES = "usingDates"
const FIX_HEIGHTS = "fixHeights"

const CORRELATION_MATRIX = "correlationMatrix"
const INVERT = "invert"

const FIRE_PARAMETER_CHANGED = "fireParameterChanged"

const VARIANCE_PROPORTION_STATISTIC = "varianceProportionStatistic"
const MATRIX_RATIO = "matrixRatio"
const COHERITABILITY = "coheritability"
const USE_EMPIRICAL_VARIANCE = "useEmpiricalVariance"

const FULL_EVALUATION = "fullEvaluation"

const MATRIX_VALIDATION = "matrixValidation"
const TRUE_PARAMETER = "trueParameter"
const INFERRED_PARAMETER = "inferredParameter"

const MODEL_EXTENSION_TRAIT_LOGGER = "modelExtensionTraitLogger"
const DIMENSIONS = "dimensions"
const ALL = "all"
const MISSING = "missing"

const INTEGRATED_FACTORS = "integratedFactorModel"
const DEFAULT_IF_NAME = "factorModel"
const LOADINGS = "loadings"
const DEFAULT_LOADINGS_ID = "L"
const PRECISION = "precision"
const FACTOR_PRECISION = "factorPrecision"

const DISTRIBUTION_LIKELIHOOD = "distributionLikelihood"
const NORMAL_DISTRIBUTION_MODEL = "normalDistributionModel"
const MEAN = "mean"
const STDEV = "stdev"

const HMC = "hamiltonianMonteCarloOperator"
const GRAD_CHECK_COUNT = "gradientCheckCount"
const JOINT_GRADIENT = "jointGradient"
const GRADIENT = "gradient"
const INTEGRATED_FACTOR_GRADIENT = "integratedFactorAnalysisLoadingsGradient"

const NSTEPS = "nSteps"
const STEP_SIZE = "stepSize"
const DRAW_VARIANCE = "drawVariance"
const GRAD_TOLERANCE = "gradientCheckTolerance"

const GAMMA_PRIOR = "gammaPrior"
const SCALE = "scale"
const SHAPE = "shape"

const NGP_OPERATOR = "normalGammaPrecisionGibbsOperator"
const NORMAL_EXTENSION = "normalExtension"
const TREE_TRAIT_NAME = "treeTraitName"

const DEFAULT_FACTOR_NAME = "factors"

const LOADINGS_GIBBS_OP = "loadingsGibbsOperator"
const RANDOM_SCAN = "randomScan"
const NEW_MODE = "newMode"
const CONSTRAINT = "constraint"
const SPARSITY = "sparsity"
const NONE = "none"
const NORMAL_PRIOR = "normalPrior"

const LATENT_FACTOR_MODEL = "latentFactorModel"
const FACTOR_NUMBER = "factorNumber"
const SCALE_DATA = "scaleData"
const FACTORS = "factors"
const DATA_FROM_TIPS = "dataAndMissingFromTreeTips"
const ROW_PRECISION = "rowPrecision"
const COL_PRECISION = "columnPrecision"
const LAT_FAC_PREC_OP = "latentFactorModelPrecisionOperator"
const FACTOR_TREE_GIBBS_OP = "factorTreeGibbsOperator"

const MULTIVARIATE_TRAIT_LIKELIHOOD = "multivariateTraitLikelihood"

const TRAIT_LOGGER = "traitLogger"
const TAXON_EXPLICIT = "taxonNameExplicit"
const NODES = "nodes"
const EXTERNAL = "external"

const PRODUCT_PARAMETER = "productParameter"
const TRANSFORMED_PARAMETER = "transformedParameter"
const TYPE = "type"
const POWER = "power"
const POWER_TRANSFORM = "powerTransform"
const TRANSFORM = "transform"

const BAYESIAN_BRIDGE = "bayesianBridge"
const GLOBAL_SCALE = "globalScale"
const EXPONENT = "exponent"
const LOCAL_SCALE = "localScale"

const MATRIX_SHRINKAGE_LIKELIHOOD = "matrixShrinkageLikelihood"
const ROW_PRIORS = "rowPriors"

const SCALE_OPERATOR = "scaleOperator"
const SCALE_FACTOR = "scaleFactor"

const BAYESIAN_BRIDGE_GIBBS_OP = "bayesianBridgeGibbsOperator"

const FACTOR_VALIDATION = "factorValidation"
const BIAS = "bias"
const SQUARED_ERROR = "squaredError"

const TRAIT_DATA_PROVIDER = "traitDataProvider"

const MULTIPLICATIVE_GAMMA_GIBBS = "multiplicativeGammaGibbsProvider"
const ROW_MULTIPLIERS = "rowMultipliers"
const ROW = "row"

const CACHED_MATRIX_INVERSE = "cachedMatrixInverse"
const COMPOUND_SYMMETRIC_MATRIX = "compoundSymmetricMatrix"
const AS_CORRELATION = "asCorrelation"
const IS_CHOLESKY = "isCholesky"
const DIAGONAL = "diagonal"
const OFF_DIAGONAL = "offDiagonal"

const LJK_CORRELATION_PRIOR = "LKJCorrelationPrior"
const SHAPE_PARAMETER = "shapeParameter"
const CHOLESKY = "cholesky"

const LOG_NORMAL_PRIOR = "logNormalPrior"
const MEAN_IN_REAL_SPACE = "meanInRealSpace"
const OFFSET = "offset"

const NORMAL_ORTHOGONAL_SUBSPACE = "normalOrthogonalSubspaceDistribution"

const NORMAL_MATRIX_NORM = "normalMatrixNormLikelihood"
const GLOBAL_PRECISION = "globalPrecision"
const MATRIX = "matrix"

const SCALED_MATRIX = "scaledMatrixParameter"
const MULTIVARIATE_GAMMA_LIKELIHOOD = "multivariateGammaLikelihood"

const MULTIPLICATIVE_PARAMETER = "multiplicativeParameter"

const GEODESIC_HMC = "geodesicHamiltonianMonteCarloOperator"
const NORMALIZED_LOADINGS_GRADIENT = "normalizedIntegratedLoadingsGradient"
const SCALE_LOADINGS_GRADIENT = "scaleIntegratedLoadingsGradient"

const MASK = "mask"
const MASKED_PARAMETER = "maskedParameter"

const INDICES = "indices"

const JOINT_OPERATOR = "jointOperator"
const INTEGRATED_FACTOR_GIBBS = "integratedFactorsGibbsOperator"
const SCALED_MATRIX_GRADIENT = "scaledMatrixGradient"
const COMPONENT = "component"
const SAMPLED_LOADINGS_GRADIENT = "sampledLoadingsGradient"

const PRECONDITIONING = "preconditioning"
const PRECONDITIONING_UPDATE_FREQUENCY = "preconditioningUpdateFrequency"
const PRECONDITIONING_MAX_UPDATE = "preconditioningMaxUpdate"
const PRECONDITIONING_DELAY = "preconditioningDelay"
const PRECONDITIONING_MEMORY = "preconditioningMemory"

const ADAPTIVE_DIAGONAL = "adaptiveDiagonal"
const MISSING_INDICATOR = "missingIndicator"

const LOCADINGS_SCALE_GIBBS = "loadingsScaleGibbsOperator"
const INDEPENDENT_NORMAL_DISTRIBUTION = "independentNormalDistributionModel"

const REJECTION_OPERATOR = "rejectionOperator"
const CONDITION = "condition"

end
