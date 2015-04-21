// -------------------------------------------------------------------------------------------------
// Test a list of valid absvars
// -------------------------------------------------------------------------------------------------

sysuse auto, clear
set more off
cls

local absvar1		turn
local absvar2		i.turn
local absvar3		tur
local absvar4		i.turn trunk
local absvar5		trunk i.turn
local absvar6		turn#trunk
local absvar7		i.turn#trunk
local absvar8		trunk#i.trunk
local absvar9		i.turn#i.trunk
local absvar10		i.turn#foreign#trunk
local absvar11		turn for trunk disp mpg
local absvar12		turn trunk#c.gear
local absvar13		turn trunk##c.gear
local absvar14		turn foreign trunk#c.(gear)
local absvar15		turn##c.(gear weight)
local absvar16		turn trunk#c.(gear weight length)
local absvar17		turn trunk#c.(gear weight length) foreign
local absvar18		turn c.(gear weight length)#trunk
local absvar19		turn c.(gear weight length)#i.trunk
local absvar20		turn (c.gear c.weight c.length)#trunk
local absvar21		turn trunk#c.gear foreign##c.(weight length)
local absvar22		FE1=turn foreign FE3=trunk
local absvar23		turn FE=foreign#c.gear
local absvar24		turn FE=foreign#c.(gear length)
local absvar25		turn foreign, savefe
local absvar26		
local absvar27		
local absvar28		
local absvar29		
local absvar30		


forval i = 1/30 {
	local absvar `absvar`i''
	cap noi ParseAbsvars `absvar'
	if _rc {
		di as error "error parsing absvar <`absvar'>"
		exit 2000
	}
}
