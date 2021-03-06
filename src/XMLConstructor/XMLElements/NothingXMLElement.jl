mutable struct NothingXMLElement <: MyXMLElement
    el::Nothing
    NothingXMLElement() = new(nothing)
end

function make_xml(::NothingXMLElement)
    # do nothing
end
