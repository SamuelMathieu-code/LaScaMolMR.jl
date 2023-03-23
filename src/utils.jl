

"""
Inverse logit function
"""
inv_logit(x) = exp(x)/(1+exp(x))


"""
Exportable Enum
"""
macro exported_enum(name, args...)
    esc(quote
        @enum($name, $(args...))
        export $name
        $([:(export $arg) for arg in args]...)
        end)
end