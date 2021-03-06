mutable struct TraitLoggerXMLElement <: MyXMLElement
    el::XMLOrNothing
    treexml::TreeModelXMLElement
    traitxml::TraitLikelihoodXMLElement

    function TraitLoggerXMLElement(treexml::TreeModelXMLElement,
                                traitxml::TraitLikelihoodXMLElement)
        return new(nothing, treexml, traitxml)
    end

end


function make_xml(xml::TraitLoggerXMLElement)
    tree_el = make_xml(xml.treexml)
    trait_el = make_xml(xml.traitxml)

    el = new_element(bn.TRAIT_LOGGER)
    trait_name = attribute(trait_el, bn.TRAIT_NAME)
    set_attributes(el, [bn.ID => bn.TRAIT_LOGGER,
                    bn.TRAIT_NAME => trait_name,
                    bn.TAXON_EXPLICIT => bn.TRUE,
                    bn.NODES => bn.EXTERNAL])


    add_ref_el(el, tree_el)
    add_ref_el(el, trait_el)
    xml.el = el
    return el
end

function get_loggables(xml::TraitLoggerXMLElement)
    make_xml(xml)
    return [xml.el]
end
