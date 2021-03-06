mutable struct CompoundXMLElement <: MyXMLElement
    els::Vector{MyXMLElement}
end

function make_xml(cp::CompoundXMLElement)
    for xml in cp.els
        make_xml(xml)
    end
end

