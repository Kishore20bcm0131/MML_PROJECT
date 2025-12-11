<!doctype html>
<html>
<head>
	<meta charset="utf-8">
	<title>Upload</title>
</head>

<body>
	<?php
	if ( ( $_FILES[ 'my_file' ][ 'name' ] != "" ) ) {
		$target_dir = "uploads/";
		$file = $_FILES[ 'my_file' ][ 'name' ];
		$path = pathinfo( $file );
		$filename = $path[ 'filename' ];
		$ext = $path[ 'extension' ];
		$temp_name = $_FILES[ 'my_file' ][ 'tmp_name' ];
		$path_filename_ext = $target_dir . $filename . "." . $ext;

		// Check if file already exists
		if ( file_exists( $path_filename_ext ) ) {
			echo "Sorry, this file name already exists.";
		} else {
			move_uploaded_file( $temp_name, $path_filename_ext );
			echo "Congratulations! File Uploaded Successfully.";
		}
	}
	?>
</body>
</html>