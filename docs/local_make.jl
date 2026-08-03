using Revise
using WannierExcitonModel
using LiveServer

# default build dir
const BUILDDIR = joinpath(@__DIR__, "build")

# 全局句柄，保证只启一个 server
const SERVER = Ref{Any}()

"""
增量生成文档并热重载浏览器。
"""
function local_make()
	# 1. 生成 html
	@info "Running makedocs..."
	include(joinpath(@__DIR__, "make.jl"))
	# 2. 第一次启动 LiveServer; 后续无需重新启动 LiveServer.
	if !isdefined(SERVER, 1) || istaskdone(SERVER, 1)
		SERVER[] = Threads.@spawn LiveServer.serve(dir = BUILDDIR, launch_browser = true)
		@info "LiveServer started."
	else
		@info "Docs rebuilt; browser will auto-refresh."
	end
	nothing
end

# 不要重复运行local_make.jl，只需要重复运行该函数或者make.jl
local_make()

