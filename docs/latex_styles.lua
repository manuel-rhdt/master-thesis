if FORMAT:match 'latex' then
  function Span(el)
    if el.classes:includes("check", 1) then
      return {
        pandoc.RawInline('latex', '\\pandoccheck{'),
        el,
        pandoc.RawInline('latex', '}')
      }
    else
      return el
    end
  end

  function Div(el)
    if el.identifier == "refs" then
      return {
        pandoc.RawBlock('latex', '{'),
        pandoc.RawBlock('latex', '\\setlength{\\parindent}{0cm}'),
        el,
        pandoc.RawBlock('latex', '}')
      }
    else
      return el
    end
  end
end