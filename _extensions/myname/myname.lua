local function makeLink(url, label)
  return string.format('<a href="%s" target="_blank">%s</a>', url, label)
end

function pathlink(webpath)
  local var = quarto.doc.input_file
  local basename = string.match(var, "[^/\\]+$")
  local namer = basename
  if next(webpath) ~= nil then
    namer = pandoc.utils.stringify(webpath[1]) .. '/' .. namer
  end
  local url = "https://github.com/chrissl789/modern-regression/raw/main/" .. namer
--  quarto.doc.includeText("before-body", makeLink(url, 'link'))
  return pandoc.RawBlock('html', makeLink(url, 'link'))
end
