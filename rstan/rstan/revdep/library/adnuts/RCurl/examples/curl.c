#include <curl/curl.h>

char *DefaultURL = "http://www.omegahat.net/index.html";


int
main(int argc, char *argv[])
{

	CURL *h;
	char **url;
	CURLcode status;

	*url = DefaultURL;

	h = curl_easy_init();
	status = curl_easy_setopt(h, CURLOPT_URL, NULL);
	if(status) {
		fprintf(stderr, "Expected error %d", status);fflush(stderr);
	}
	curl_easy_setopt(h, CURLOPT_URL, *url);
	curl_easy_perform(h);

	return(0);
}
