mutable struct TimerXMLElement <: MyXMLElement
    el::XMLOrNothing
    mcmc_el::MCMCXMLElement
    filename::String

    TimerXMLElement(mcmc_el::MCMCXMLElement) = new(nothing, mcmc_el, "")
    TimerXMLElement(mcmc_el::MCMCXMLElement, filename::String) = new(nothing, mcmc_el, filename)
end

function make_xml(t_el::TimerXMLElement)
    el = new_element(bn.REPORT)

    if length(t_el.filename) > 0
        set_attribute(el, bn.FILENAME, t_el.filename)
    end

    prop_el = new_child(el, bn.PROPERTY)
    set_attribute(prop_el, bn.NAME, bn.TIMER)
    add_ref_el(prop_el, t_el.mcmc_el.el)

    t_el.el = el
    return el
end

function set_filename(t::TimerXMLElement, fn::String)
    t.filename = fn
end
