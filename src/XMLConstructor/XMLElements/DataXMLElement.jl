mutable struct DataXMLElement <: MyXMLElement
    el::XMLOrNothing
    data_mats::Vector{<:AbstractMatrix{Float64}}
    trait_names::Vector{S} where S <: AbstractString
    taxa::Vector{<:AbstractString}
    node_times::Vector{Float64}
    use_times::Bool

    function DataXMLElement(mat::AbstractMatrix{Float64},
                taxa::Vector{S},
                newick::String;
                trait_name::AbstractString = bn.DEFAULT_TRAIT_NAME
                ) where S <: AbstractString

        node_times = get_node_times(taxa, newick)
        return new(nothing, [mat], [trait_name], taxa, node_times, false)
    end

    function DataXMLElement(mats::Vector{<:AbstractMatrix{Float64}},
                trait_names::Vector{T} where T <: AbstractString,
                taxa::Vector{S},
                newick::String
                ) where S <: AbstractString

        node_times = get_node_times(taxa, newick)

        return new(nothing, mats, trait_names, taxa, node_times, false)
    end

end

function DataXMLElement(dm::DataModel, newick::String)
    n_partitions = length(dm.partitions)
    n = size(dm.data, 1)
    data = [@view dm.data[1:n, dm.partitions[i]] for i = 1:n_partitions]
    return DataXMLElement(data, dm.trait_names, dm.taxa, newick)
end




function make_xml(dl::DataXMLElement)
    dl.el = make_taxa_el(dl)
end

function make_taxa_el(dxl::DataXMLElement)

    data_mats = dxl.data_mats
    el = new_element(bn.TAXA)
    n, p = size(data_mats[1])
    for i = 1:length(data_mats)
        @assert size(data_mats[i], 1) == n
    end
    for i = 1:n
        add_taxon(dxl, el, i)
    end
    return el
end

function add_taxon(dxl::DataXMLElement, pel::XMLElement, ind::Int)

    taxa = dxl.taxa
    data_mats = dxl.data_mats
    trait_names = dxl.trait_names
    node_times = dxl.node_times


    el = new_element(bn.TAXON)
    set_attribute(el, bn.ID, taxa[ind])
    if dxl.use_times
        date_el = new_child(el, bn.DATE)
        set_attributes(date_el, [(bn.VALUE, string(node_times[ind])), (bn.DIRECTION, bn.FORWARDS)])
    end
    for i = 1:length(data_mats)
        attr_el = new_element(bn.ATTR)
        set_attribute(attr_el, bn.NAME, trait_names[i])
        s_data = format_data_line(data_mats[i], ind)
        set_content(attr_el, s_data)
        add_child(el, attr_el)
    end
    add_child(pel, el)
    return el
end

function format_data_line(data::AbstractMatrix{Float64}, ind::Int)
    p = size(data, 2)
    s_data = Vector{String}(undef, p)
    for j = 1:p
        val = data[ind, j]
        if isnan(val)
            s_data[j] = bn.MISSING_VAL
        else
            s_data[j] = string(val)
        end
    end
    return join(s_data, ' ')
end

function get_node_times(taxa::Vector{String}, newick::String)
    tree = PhyloNetworks.readTopology(newick)
    root_dists = TreeUtils.leaf_distances(tree)
    tree_taxa = [tree.leaf[i].name for i = 1:tree.numTaxa]
    node_times = zeros(tree.numTaxa)
    for i = 1:length(taxa)
        taxon = taxa[i]
        ind = findfirst(x -> x == taxon, tree_taxa)
        node_times[i] = root_dists[ind]
    end
    return node_times
end

function add_trait!(dxl::DataXMLElement, data::Matrix{Float64}, trait_name::String)
    push!(dxl.data_mats, data)
    push!(dxl.trait_names, trait_name)
end

function set_dates(xml::DataXMLElement, dates::AbstractVector{Float64})
    if length(dates) != length(xml.node_times)
        DimensionMismatch("The number of dates must match the number of taxa.")
    end

    xml.node_times .= dates
end