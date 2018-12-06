using Documenter, NonparEconometricsTool

makedocs(;
    modules=[NonparEconometricsTool],
    format=:html,
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/JasmineHao/NonparEconometricsTool.jl/blob/{commit}{path}#L{line}",
    sitename="NonparEconometricsTool.jl",
    authors="Jasmine Hao, University of British Columbia",
    assets=[],
)

deploydocs(;
    repo="github.com/JasmineHao/NonparEconometricsTool.jl",
    target="build",
    julia="1.0",
    deps=nothing,
    make=nothing,
)
