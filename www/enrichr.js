// http://amp.pharm.mssm.edu/Enrichr/help#faq
Shiny.addCustomMessageHandler("enrich",
	function(options) {
	    if (typeof options.list === 'undefined') {
	        alert('No genes defined.');
	    }

	    var description  = options.description || "",
	    	popup = options.popup || false,
	    	form = document.createElement('form'),
	    	listField = document.createElement('input'),
	    	descField = document.createElement('input');
	  
	    form.setAttribute('method', 'post');
	    form.setAttribute('action', 'http://amp.pharm.mssm.edu/Enrichr/enrich');
	    if (popup) {
	        form.setAttribute('target', '_blank');
	    }
	    form.setAttribute('enctype', 'multipart/form-data');

	    listField.setAttribute('type', 'hidden');
	    listField.setAttribute('name', 'list');
	    listField.setAttribute('value', options.list);
	    form.appendChild(listField);

	    descField.setAttribute('type', 'hidden');
	    descField.setAttribute('name', 'description');
	    descField.setAttribute('value', description);
	    form.appendChild(descField);

	    document.body.appendChild(form);
	    form.submit();
	    document.body.removeChild(form);
	}
);