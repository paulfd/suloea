
function (event) {
	if (event.type == 'change' && event.symbol == 'retuning') {
		if (event.value == 1) {
			event.icon.find("[mod-role=retuning]").each(function () { $(this).removeClass("off"); });
			event.icon.find("[mod-role=retuning]").each(function () { $(this).addClass("on"); });
		} else {
			event.icon.find("[mod-role=retuning]").each(function () { $(this).removeClass("on"); });
			event.icon.find("[mod-role=retuning]").each(function () { $(this).addClass("off"); });
		}
	}
}