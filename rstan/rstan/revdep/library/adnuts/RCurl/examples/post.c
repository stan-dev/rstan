#include <curl/curl.h>


int
main(int argc, char *argv[])
{

	CURL *h;
	char url[] = "http://www.speakeasy.org/~cgires/perl_form.cgi";
	CURLcode status;

	struct curl_httppost *post = NULL, *last = NULL;

	h = curl_easy_init();

	curl_easy_setopt(h, CURLOPT_URL, url);


	curl_formadd(&post, &last, 
		     CURLFORM_PTRNAME, "some_text",
		     CURLFORM_PTRCONTENTS, "Duncan",
		     CURLFORM_END);

	curl_formadd(&post, &last, 
		     CURLFORM_PTRNAME, "choice",
		     CURLFORM_PTRCONTENTS, "Ho",
		     CURLFORM_END);

	curl_formadd(&post, &last, 
		     CURLFORM_PTRNAME, "radbut",
		     CURLFORM_PTRCONTENTS, "eep",
		     CURLFORM_END);

	curl_formadd(&post, &last, 
		     CURLFORM_PTRNAME, "box",
		     CURLFORM_PTRCONTENTS, "box1",
		     CURLFORM_END);

        curl_easy_setopt(h, CURLOPT_HTTPPOST, post);
	curl_easy_perform(h);

	return(0);
}
