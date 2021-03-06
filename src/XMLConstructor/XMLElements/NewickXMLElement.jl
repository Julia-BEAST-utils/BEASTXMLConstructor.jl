mutable struct NewickXMLElement <: MyXMLElement
    el::XMLOrNothing
    newick::String
    fix_tree::Bool

    NewickXMLElement(newick::String) = new(nothing, newick, true)
end

function make_xml(nl::NewickXMLElement)

    el = new_element(bn.NEWICK)
    set_attribute(el, bn.ID, bn.DEFAULT_TREE_NAME)
    if nl.fix_tree
        set_attribute(el, bn.USING_HEIGHTS, bn.TRUE)
        set_attribute(el, bn.USING_DATES, bn.FALSE)
    end
    set_content(el, "$(nl.newick);")
    nl.el = el
    return el
end
