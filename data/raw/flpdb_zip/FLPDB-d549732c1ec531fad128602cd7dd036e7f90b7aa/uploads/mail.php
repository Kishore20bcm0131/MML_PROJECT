<?php
$filenameee = $_FILES[ 'upload' ][ 'name' ];
$fileName = $_FILES[ 'upload' ][ 'tmp_name' ];
$firstname = $_POST[ 'firstname' ];
$lastname = $_POST[ 'lastname' ];
$email = $_POST[ 'email' ];
$orcid = $_POST[ 'orcid' ];
$institution = $_POST[ 'institution' ];
$moreinfo = $_POST[ 'info' ];

$message = "First Name = " . $firstname . "\r\n Last Name = " . $lastname . "\r\n  Email = " . $email . "\r\n  ORCID = " . $orcid . "\r\n  Institution = " . $institution . "\r\n More info =" . $moreinfo;

$subject = "FLP File Submitted";
$fromname = "FLP Database";
$fromemail = 'nobody@clarkson.edu';

$mailto = 'jye@clarkson.edu';

if ( !empty( $_FILES[ "upload" ][ "name" ] ) ) {

	// File path config
	$targetDir = "uploads/";
	$fileNamee = basename( $_FILES[ "upload" ][ "name" ] );
	$targetFilePath = $targetDir . $fileNamee;
	$fileType = pathinfo( $targetFilePath, PATHINFO_EXTENSION );

	// Allow certain file formats
	$allowTypes = array( 'cif', 'xyz' );
	if ( in_array( $fileType, $allowTypes ) ) {
		// Upload file to the server
		if ( move_uploaded_file( $_FILES[ "upload" ][ "tmp_name" ], $targetFilePath ) ) {
			$uploadedFile = $targetFilePath;
			$uploadStatus = 1;

		} else {
			$uploadStatus = 0;
			echo "Sorry, there was an error uploading your file.";
		}
	} else {
		$uploadStatus = 0;
		echo 'Sorry, only CIF, & XYZ files are allowed to upload.';
	}
}
if ( $uploadStatus == 1 ) {

	$content = file_get_contents( $fileName );
	$content = chunk_split( base64_encode( $content ) );
	// a random hash will be necessary to send mixed content
	$separator = md5( time() );
	// carriage return type (RFC)
	$eol = "\r\n";
	// main header (multipart mandatory)
	$headers = "From: " . $fromname . " <" . $fromemail . ">" . $eol;
	$headers .= "MIME-Version: 1.0" . $eol;
	$headers .= "Content-Type: multipart/mixed; boundary=\"" . $separator . "\"" . $eol;
	$headers .= "Content-Transfer-Encoding: 7bit" . $eol;
	$headers .= "This is a MIME encoded message." . $eol;
	// message
	$body = "--" . $separator . $eol;
	$body .= "Content-Type: text/plain; charset=\"iso-8859-1\"" . $eol;
	$body .= "Content-Transfer-Encoding: 8bit" . $eol;
	$body .= $message . $eol;
	// attachment
	$body .= "--" . $separator . $eol;
	$body .= "Content-Type: application/octet-stream; name=\"" . $filenameee . "\"" . $eol;
	$body .= "Content-Transfer-Encoding: base64" . $eol;
	$body .= "Content-Disposition: attachment" . $eol;
	$body .= $content . $eol;
	$body .= "--" . $separator . "--";
	//SEND Mail
	if ( mail( $mailto, $subject, $body, $headers ) ) {
		echo "Your file has been submitted successfully!";
		// Delete attachment file from the server
		@unlink( $uploadedFile );


	} else {
		echo "Your file submission failed, please try again.";
		print_r( error_get_last() );
		@unlink( $uploadedFile );
	}
}