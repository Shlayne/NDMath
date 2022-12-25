require('vstudio')

local p = premake
local m = p.vstudio.vc2010

p.api.register {
	name = "usestdpreproc",
	scope = "config",
	kind = "boolean"
}

p.api.register {
	name = "enablemodules",
	scope = "config",
	kind = "boolean"
}

p.override(m.elements, "clCompile", function(base, cfg)
	local calls = base(cfg)

	table.replace(calls, m.floatingPointModel, function(cfg)
		if cfg.floatingpoint and cfg.floatingpoint ~= "Default" then
			if cfg.floatingpoint ~= "None" then
				m.element("FloatingPointModel", nil, cfg.floatingpoint)
			else
				p.w("<FloatingPointModel />") -- Having "None" as the value causes visual studio's build to crash.
			end
		end
	end)

	table.insertafter(calls, m.languageStandardC, function(cfg)
		if _ACTION >= "vs2019" and cfg.usestdpreproc ~= nil then
			m.element("UseStandardPreprocessor", nil, iif(cfg.usestdpreproc, "true", "false"))
		end
	end)
	
	table.insertafter(calls, m.scanSourceForModuleDependencies, function(cfg)
		if _ACTION >= "vs2022" and cfg.enablemodules ~= nil then
			m.element("EnableModules", nil, iif(cfg.enablemodules, "true", "false"))
		end
	end)

	return calls
end)
