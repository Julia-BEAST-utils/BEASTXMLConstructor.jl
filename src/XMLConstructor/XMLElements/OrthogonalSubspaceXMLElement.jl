mutable struct OrthogonalSubspaceXMLElement <: MyXMLElement
    el::XMLOrNothing
    loadings::MatrixParameter
    prior::MatrixShrinkageLikelihoods

    function OrthogonalSubspaceXMLElement(loadings::MatrixParameter,
                                          prior::MyXMLElement)
        return new(nothing, loadings, prior)
    end
end



