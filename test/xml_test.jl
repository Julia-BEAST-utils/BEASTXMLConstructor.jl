using BEASTXMLConstructor
using BeastUtils.TreeUtils, PhyloNetworks

xc = BEASTXMLConstructor

using DataFrames, CSV, Test
### Making sure https://github.com/suchard-group/incomplete_measurements/blob/master/scripts/xml_setup.jl works

function make_mbd(data_path::String, newick_path::String, xml_path::String,
                filename::String; dates_path::String = "")

    df = DataFrame(CSV.File(data_path))

    use_dates = false
    if length(dates_path) > 0
        dates_df = DataFrame(CSV.File(dates_path))
        @assert dates_df[!, :taxon] == df[!, :taxon]
        use_dates = true
    end

    newick = read(newick_path, String)

    taxa, data = BEASTXMLConstructor.df_to_matrix(df)


    bx = BEASTXMLConstructor.make_residual_xml(data, taxa, newick, chain_length = 100_000)
    if use_dates
        BEASTXMLConstructor.use_dates!(bx)
        xc.set_data_dates(bx, dates_df[!, :date])
    end
    xc.set_screen_logEvery(bx, 100)
    xc.set_file_logEvery(bx, 10)
    xc.set_filename(bx, filename)
    BEASTXMLConstructor.add_MBD_loggables!(bx)

    BEASTXMLConstructor.save_xml(xml_path, bx)
end


n = 10
p = 4
taxa = ["taxon$i" for i = 1:n]
data = randn(n, p)
dates = rand(n)

df = DataFrame()
df.taxon = taxa
for i = 1:p
    df[!, Symbol("trait$i")] = data[:, i]
end

data_path = "data.csv"
dates_path = "dates.csv"
newick_path = "newick.txt"
xml_path = "BEASTXMLConstructor.xml"
dates_xml_path = "xml_dates.xml"
filename = "test"


CSV.write(data_path, df)
CSV.write(dates_path, DataFrame(taxon = taxa, date = dates))

newick = writeTopology(TreeUtils.rtree(taxa))
write(newick_path, newick)


make_mbd(data_path, newick_path, xml_path, filename)

@test isfile(xml_path)

make_mbd(data_path, newick_path, dates_xml_path, filename, dates_path = dates_path)

@test isfile(dates_xml_path)



### Testing Factor xml
k = 3

# HMC, no shrinkage
bx = BEASTXMLConstructor.make_pfa_xml(data, taxa, newick, k, useHMC = false)
BEASTXMLConstructor.save_xml("facHMC.xml", bx)
@test isfile("facHMC.xml")

# Gibbs, no shrinkage
bx = BEASTXMLConstructor.make_pfa_xml(data, taxa, newick, k, useHMC = true)
BEASTXMLConstructor.save_xml("facGibbs.xml", bx)
@test isfile("facGibbs.xml")

# # HMC, shrinkage
# bx = BEASTXMLConstructor.make_pfa_xml(data, taxa, newick, k, useHMC = true,
#             shrink_loadings = true)
# BEASTXMLConstructor.save_xml("facHMCShrink.xml", bx)
# @test isfile("facHMCShrink.xml")

# # Gibbs, shrinkage
# bx = BEASTXMLConstructor.make_pfa_xml(data, taxa, newick, k, useHMC = false,
#             shrink_loadings = true)
# BEASTXMLConstructor.save_xml("facGibbsShrink.xml", bx)
# @test isfile("facGibbsShrink.xml")

# # Rotate prior
# fn = "facGibbsRotate.xml"
# bx = BEASTXMLConstructor.make_pfa_xml(data, taxa, newick, k, useHMC = false,
#             shrink_loadings = true, rotate_prior = true)
# BEASTXMLConstructor.save_xml(fn, bx)
# @test isfile(fn)

################################################################################
## Orthogonal factor analysis
################################################################################

fn = "facOrthogonal.xml"
bx = BEASTXMLConstructor.make_orthogonal_pfa_xml(data, taxa, newick, k)
BEASTXMLConstructor.save_xml(fn, bx)
@test isfile(fn)


################################################################################
## Joint models
################################################################################

n = 10
k = 2
p_res = 3
p_fac = 5
p_diff = k + p_res

taxa = ["taxon_$i" for i = 1:n]

tree = rtree(taxa)
newick = writeTopology(tree)

rm = ResidualVarianceModel(p_res)
fm = IntegratedFactorModel(k, p_fac)

jm = JointProcessModel([rm, fm])

dm = DataModel(taxa,
               [randn(n, p_res), randn(n, p_fac)],
               ["trait.res", "trait.fac"])

bx = make_joint_xml(newick, dm, jm)
joint_path = "joint.xml"
# BEASTXMLConstructor.save_xml(joint_path, bx)
# @test isfila(joint_path)
# rm(joint_path)

################################################################################
## Integrate two xml
################################################################################

test_path = joinpath(@__DIR__, "data", "sequence.xml")

bx_seq = BEASTXMLElement(test_path)

taxa = BEASTXMLConstructor.find_element(bx_seq, BEASTXMLConstructor.EmptyDataXMLElement).taxa

k = 2
p = 10
n = length(taxa)
data = randn(n, p)
newick = writeTopology(rtree(taxa))

bx = BEASTXMLConstructor.make_pfa_xml(data, taxa, newick, k, useHMC = false)
BEASTXMLConstructor.save_xml("facGibbs.xml", bx)
@test isfile("facGibbs.xml")

BEASTXMLConstructor.merge_xml!(bx, bx_seq)
# xml = BEASTXMLConstructor.make_xml(bx)
BEASTXMLConstructor.save_xml("merge.xml", bx);
# print(xml)

################################################################################
## Sampled factor model with integrated factors sampler
################################################################################
k = 2
p = 10
n = length(taxa)
data = randn(n, p)
newick = writeTopology(rtree(taxa))

bx = BEASTXMLConstructor.make_sampled_pfa_xml(data, taxa, newick, k)
save_xml("facSampledHMC.xml", bx, change_filename=true)
