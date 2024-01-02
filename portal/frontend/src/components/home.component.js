import React, { Component } from "react";
import { Navigate } from "react-router-dom";

import AuthService from "../services/auth.service";
import UserService from "../services/user.service";

export default class Home extends Component {
  constructor(props) {
    super(props);

    this.state = {
      content: ""
    };
  }

  componentDidMount() {
    const currentUser = AuthService.getCurrentUser();
    if (!currentUser) this.setState({ redirect: "/login" });
    this.setState({ currentUser: currentUser, userReady: true })

    UserService.getPublicContent().then(
      response => {
        this.setState({
          content: response.data
        });
      },
      error => {
        this.setState({ redirect: "/login" });
      }
    );
  }

  render() {
    if (this.state.redirect) {
      return <Navigate to={this.state.redirect} />
    }

    return (
      <iframe 
        src={this.state.content}// URL of the content to be loaded
        width="100%"  // Width of the iframe
        height="600"  // Height of the iframe
        frameBorder="0" // Style of the border, set to 0 for no border
        allowFullScreen // Allow full screen if the content supports it
      >
        {/* Fallback content in case iframes are not supported */}
        Your browser does not support iframes.
      </iframe>
    );
  }
}